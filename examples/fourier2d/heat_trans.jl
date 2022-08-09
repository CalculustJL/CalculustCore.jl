#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, LinearAlgebra
using Plots, Test

nx = 32
ny = 32
ν = 1e-2
p = nothing

""" space discr """
space  = FourierSpace(nx, ny)
tspace = transform(space)
discr  = Collocation()

x, y = points(space)
(k,) = points(tspace)
F    = transformOp(space)

α = 5
β = 3
u0 = @. sin(α*x) * sin(β*y)
û0 = F * u0

Â = diffusionOp(ν, tspace, discr)
F̂ = SciMLOperators.NullOperator(tspace)

odefunc = cache_operator(Â + F̂, k)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Tsit5()
prob = ODEProblem(odefunc, û0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave, abstol=1e-8, reltol=1e-8)

""" analysis """
pred = [F,] .\ sol.u
pred = hcat(pred...)

utrue(t) = @. u0 * (exp(-ν*(α^2+β^2)*t))
ut = utrue(sol.t[1])
for i=2:length(sol.t)
    utt = utrue(sol.t[i])
    global ut = hcat(ut, utt)
end

err = norm(pred .- ut, Inf)
println("frobenius norm of error across time", err)
@test err < 1e-7
#
