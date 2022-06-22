#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces
using OrdinaryDiffEq, LinearSolve
using Plots

N = 1024
ν = 0e-2
p = ()

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)
u0 = @. sin(1x)

A = diffusionOp(ν, space, discr)

# verify update behaviour
function burgers!(L, u, p, t)
    L.diag .= u
    L
end

v = @. x*0 + 1
f = @. x*0

C = advectionOp((v,), space, discr)
F = AffineOperator(C, f)

#C = advectionOp((u0,), space, discr; vel_update_func=burgers!)
#F = C

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
tspan = (0.0, 10.0)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg)
@show sol.retcode

nothing
#
