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

N = 128
ν = 1e-2
p = ()

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
ftr  = transformOp(space)

α = 2
uic(x) = @. sin(α*x)
function utrue(t,x)
    cos(t) * uic(x)
end

A = diffusionOp(ν, space, discr)

f = @. x*0 + .1
function forcing!(f, u, p, t)
    ui = -sin(t)*uic(x)
    ud = -ν*α*α*uic(x)*cos(t)
    f .= ui - ud
    f
end
F = forcingOp(f, space, discr; f_update_func=forcing!)

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
u0 = uic(x)
tspan = (0.0, 2π)
tsave = range(tspan...; length=10)
odealg = Tsit5()
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave, reltol=1e-8, abstol=1e-8)

""" analysis """
pred = Array(sol)

ut = utrue(sol.t[1],x)
for i=2:length(sol.t)
    utt = utrue(sol.t[i],x)
    global ut = hcat(ut, utt)
end

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=false)
end
display(plt)

err = norm(pred .- ut, Inf)
@test err < 1e-8
#
