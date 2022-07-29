#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearAlgebra, Plots, Test

N = 128
ν = 1e-2
p = nothing

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)

α = 2
uic(x) = @. sin(α*x)
u0 = uic(x)
function utrue(t,x)
    cos(t) * uic(x)
end

A = diffusionOp(ν, space, discr)
function forcing!(f, u, p, t)
    ui = -sin(t)*uic(x)
    ud = -ν*α*α*uic(x)*cos(t)
    f .= ui - ud
    f
end
F = forcingOp(zero(u0), space, discr; f_update_func=forcing!)
ddt = cache_operator(A+F, u0)

function dudt(du, u, p, t)
    ddt(du, u, p, t)
end

""" time discr """
tspan = (0.0, 2π)
tsave = range(tspan...; length=10)

odeprob = ODEProblem(dudt, u0, tspan, p)
@time sol = solve(odeprob, Tsit5(), saveat=tsave, reltol=1e-10, abstol=1e-10)

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
display(err)
@test err < 1e-7
#
