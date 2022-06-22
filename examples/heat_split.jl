#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces
using OrdinaryDiffEq, LinearSolve
using Plots

N = 1024
ν = 1e-2
p = ()

""" space discr """
domain = FourierDomain(1)
space  = FourierSpace(N; domain=domain)
discr  = Collocation()

(x,) = points(space)

D = diffusionOp(ν, space, discr)
D = cache_operator(D, x)

F = SciMLOperators.NullOperator{length(space)}()

""" IC """
u0 = @. sin(10x)

""" time discr """
tspan = (0.0, 10.0)
prob = SplitODEProblem(D, F, u0, tspan, p)

odealg = Rodas5(autodiff=false)
@time sol  = solve(prob, odealg)
@show sol.retcode

nothing
#
