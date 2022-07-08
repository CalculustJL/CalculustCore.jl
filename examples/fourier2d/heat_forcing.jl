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

α = 2
β = 3
uic(x, y) = @. sin(α*x)*sin(β*y)
utrue(t,x, y) = cos(t) * uic(x, y)

A = diffusionOp(ν, space, discr)

function forcing!(f, u, p, t)
    ui = -sin(t)*uic(x,y)
    ud = -ν*(α^2 + β^2)*uic(x,y)*cos(t)
    f .= ui - ud
    f
end
F = forcingOp(zero(x), space, discr; f_update_func=forcing!)

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
u0 = uic(x,y)
tspan = (0.0, 2π)
tsave = range(tspan...; length=10)
odealg = Tsit5()

prob = SplitODEProblem(A, F, u0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave, reltol=1e-10, abstol=1e-10)

""" analysis """
pred = Array(sol)

ut = utrue(sol.t[1],x,y)
for i=2:length(sol.t)
    utt = utrue(sol.t[i],x,y)
    global ut = hcat(ut, utt)
end

anim = animate(pred, space)
filename = joinpath(dirname(@__FILE__), "heat_forcing" * ".gif")
gif(anim, filename, fps=5)

err = norm(pred .- ut, Inf)
display(err)
@test err < 1e-7
#
