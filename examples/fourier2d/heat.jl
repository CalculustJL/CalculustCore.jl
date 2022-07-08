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
space = FourierSpace(nx, ny)
discr = Collocation()

x, y = points(space)
ftr  = transformOp(space)

""" IC """
α = 5
β = 3
u0 = @. sin(α*x) * sin(β*y)

""" operators """
A = diffusionOp(ν, space, discr)
F = SciMLOperators.NullOperator(space)

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Tsit5()
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave, reltol=1e-8)

""" analysis """
pred = Array(sol)

utrue(t) = @. u0 * (exp(-ν*(α^2+β^2)*t))
ut = utrue(sol.t[1])
for i=2:length(sol.t)
    utt = utrue(sol.t[i])
    global ut = hcat(ut, utt)
end

anim = animate(pred, space)
filename = joinpath(dirname(@__FILE__), "heat" * ".gif")
gif(anim, filename, fps=5)

err = norm(pred .- ut, Inf)
println("frobenius norm of error across time", err)
@test err < 1e-7
#
