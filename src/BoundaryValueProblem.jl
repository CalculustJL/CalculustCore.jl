#
abstract type AbstractBoundaryValueProblem <: SciMLBase.DEProblem end
abstract type AbstractBoundaryValueAlgorithm <: SciMLBase.DEAlgorithm end

"""
lhs(u) = b

move boundary data to right hand side, and solve BVP A \ b
using LinearSolve.jl
"""
struct BoundaryValuePDEProblem{Tu,Tbc,Tlhsop,Trhs,Tspace} <: AbstractBoundaryValueProblem
    u::Tu
    bc::Tbc
    lhsOp::Tlhsop
    rhs::Trhs
    space::Tspace
end

struct LinearBVPDEAlg{Tl} <: AbstractBoundaryValueAlgorithm
    linalg::Tl
end

# integrate NonlinearSolve.jl with LinearSolve.jl first
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
    rhs = b - applyBC()
end

function SciMLBase.sovle(prob::BoundaryValuePDEProblem, alg::AbstractBoundaryValueProblem)

    b = makeRHS(prob)

    linprob = LinearProblem(lhsOp, b; u0=u)
    linsol = solve(linprob, linalg)

    linsol.u
end
#
