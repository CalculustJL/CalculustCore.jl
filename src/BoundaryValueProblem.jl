#
abstract type AbstractBoundaryValueProblem <: SciMLBase.DEProblem end
abstract type AbstractBoundaryValueAlgorithm <: SciMLBase.DEAlgorithm end

struct BoundaryValuePDEProblem{Tu,Tbc,Top,Tf,Tsp} <: AbstractBoundaryValueProblem
    """Neumann Operator"""
    op::Top
    """Right-hand-side function"""
    f::Tf
    """Initial guess"""
    u::Tu
    """Boundary condition dictionary"""
    bc::Tbc
    """Function space"""
    space::Tsp
end

function BoundaryValuePDEProblem(op, f, bcs, space; u=nothing)
    bc = BoundaryCondition(bcs, space)

    if u isa Nothing
        x = get_grid(space)[1]
        u = zero(x)
    end

    BoundaryValuePDEProblem(op, f, u, bc, space)
end

struct LinearBVPDEAlg{Tl} <: AbstractBoundaryValueAlgorithm
    linalg::Tl
end

#TODO integrate NonlinearSolve.jl with LinearSolve.jl first
struct NonlinearBVPDEAlg{Tnl} <: AbstractBoundaryValueAlgorithm
    nlalg::Tnl
end

"""
Solve boundary value problem
A * u = f

when A is a linear operator, we get

M * A * uh = M * ( B * f - A * ub)

where u = uh + ub (Dirichlet data)

learn to add neumann, robin data, and solve BVP

plug in to SciMLBase.BVProblem
"""
function makeRHS(f, bc::BoundaryCondition)

    M = MassOp(space)
    b = M * f

    dirichlet = zero(b)
    neumann   = zero(b)
    robin     = zero(b)

    grid = get_grid(space)

    for i=1:num_boundaries(domain)
        tag = boundary_tag(domain, i)
        bc  = bcs[tag]
        idx = indices[i]

        # zygote.ignore()
        bcnodes = begin
            x = get_grid(space)[1]
            M = Bool.(zero(x))
            M[idx] = true
        end
        #

        if bc isa DirichletBC
            dirichlet += bc.f(grid) * bcnodes
        elseif bc isa NeumannBC
            neumann += bc.f(grid) * bcnodes
        elseif bc isa RobinBC
            # TODO
        elseif bc isa PeriodicBC
            continue
        end
    end

    b = mask * b + dirichlet - neumann + robin
end

"""
 [A_II A_IB] * [u_I] = [b_I]
 [A_BI A_BB]   [u_B]   [b_B]

with mask M,
    u_I = M * u       # masks ∂ data
    u_B = (I - M) * u # hides interior data
"""
function SciMLBase.solve(prob::BoundaryValuePDEProblem, alg::AbstractBoundaryValueProblem)
    @unpack op, f, u, bc, space = prob

    b = makeRHS(f, bc)

    linprob = LinearProblem(op, b; u0=u)
    linsol = solve(linprob, linalg)

    linsol.u
end
#
