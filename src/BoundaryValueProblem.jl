#
abstract type AbstractBoundaryValueProblem <: SciMLBase.DEProblem end
abstract type AbstractBoundaryValueAlgorithm <: SciMLBase.DEAlgorithm end

struct BoundaryValuePDEProblem{Tu,Tbc,Top,Tf,Tsp} <: AbstractBoundaryValueProblem
    """Neumann Operator"""
    op::Top
    """Right-hand-side vector"""
    f::Tf
    """Initial guess"""
    u::Tu
    """Boundary condition object"""
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

function makeLHS(op::AbstractOperator{<:Number,D},
                 bc::AbstractBoundaryCondition{<:Number,D}) where{D}
    @unpack mask_dir, amask_dir = bc

    #TODO
    """
    how to empty dirichlet row? -
        https://github.com/vpuri3/PDEInterfaces.jl/issues/1

    is having the construct u := u_inhom + u_hom the only way?
    then would have to interpolate u_inhom into interior.
    what are the consequences?
    """

    lhs = mask_dir * op + amask_dir
end

function makeRHS(f::AbstractField{<:Number,D},
                 bc::AbstractBoundaryCondition{<:Number,D}) where{D}
    @unpack bcs, antimasks, dirichlet_mask, space = bc

    M = MassOp(space)
    b = M * f

    dirichlet = zero(b)
    neumann   = zero(b)
    robin     = zero(b)

    grid = get_grid(space)

    for i=1:num_boundaries(domain)
        tag   = boundary_tag(domain, i)
        bc    = bcs[tag]
        amask = antimasks[i]

        if bc isa DirichletBC
            dirichlet += amask * bc.f(grid...)
        elseif bc isa NeumannBC
            neumann += amask * bc.f(grid...)
        elseif bc isa RobinBC
            # TODO
        elseif bc isa PeriodicBC
            continue
        end
    end

    b = dirichlet_mask * b + dirichlet - neumann + robin
end

function SciMLBase.solve(prob::BoundaryValuePDEProblem, alg::AbstractBoundaryValueProblem)
    @unpack op, f, u, bc, space = prob
    @unpack linalg = alg

    grid = get_grid(space)
    ff = f(grid...)

    lhs = makeLHS(op, bc)
    rhs = makeRHS(ff, bc)

    linprob = LinearProblem(op, b; u0=u)
    linsol = solve(linprob, linalg)

    linsol.u
end
#
