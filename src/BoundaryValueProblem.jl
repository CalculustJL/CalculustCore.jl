#
abstract type AbstractBoundaryValueProblem <: SciMLBase.DEProblem end
abstract type AbstractBoundaryValueAlgorithm <: SciMLBase.DEAlgorithm end

struct BoundaryValuePDEProblem{Tu,Tbc,Tlhsop,Trhs,Tspace} <: AbstractBoundaryValueProblem
    u::Tu
    bc::Tbc
    op::Top         # neumann operator
    space::Tsp
    mass_matrix::Tm # RHS mass_matrix
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
function makeRHS(prob::BoundaryValuePDEProblem)
    rhs = b - applyBC(u) # dirichlet
    rhs = b + applyBC(u) # neumann
end

"""
 [A_II A_IB] * [u_I] = [b_I]
 [A_BI A_BB]   [u_B]   [b_B]

with mask M,
    u_I = M * u       # masks ∂ data
    u_B = (I - M) * u # hides interior data
"""

function SciMLBase.sovle(prob::BoundaryValuePDEProblem, alg::AbstractBoundaryValueProblem)

    b = makeRHS(prob)

    linprob = LinearProblem(lhsOp, b; u0=u)
    linsol = solve(linprob, linalg)

    linsol.u
end
#
