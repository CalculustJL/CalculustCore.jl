#
"""
Solve boundary value problem
A * u = f

where A is a linear operator
we get

M * A * uh = M * ( B * f - A * ub)

where u = uh + ub (Dirichlet data)

learn to add neumann, robin data, and solve BVP

plug in to SciMLBase.BVProblem
"""
function makeRHS(prob::BoundaryValueProblem)

    rhs = b - applyBC()
end

function SciMLBase.solve(prob::BoundaryValueProblem)
    BoundaryValueCache()
end

"""
lhs(u) = b

move boundary data to right hand side, and solve BVP A \ b
using LinearSolve.jl
"""
struct BoundaryValuePDEProblem
    u0::Tu0
    bc::Tbc

    lhs::Tlhs # op or func

    rhs_func::Trhs
    b

    space
end

struct BoundaryValuePDECache
    u0::Tu
    bc::Tbc

    lhs::Tlhs

    rhs_func
    b::Tu

    space

    linprob
end

function SciMLBase.init(prob::BoundaryValuePDEProblem, alg::SciMLBase.AbstractBoundaryValueAlgorithm)
end

#
