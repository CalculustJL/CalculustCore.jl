#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces
using OrdinaryDiffEq, Zygote, Lux
using Plots

N = 128
ν = 1.0
p = ()

""" time discr """
tspan = (0.0, 1.0)
odealg = Tsit5()

""" space discr """
domain = FourierDomain(1)
space  = FourierSpace(N; domain=domain)
discr  = Collocation()

(x,) = points(space)

D = diffusionOp(ν, space, discr)
D = cache_operator(D, x)

""" IC """
u0 = @. sin(x)

""" solve """
prob = ODEProblem(D, u0, tspan, p)
sol  = solve(prob, odealg)
