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
ν = 0e0
p = nothing

""" space discr """
space = FourierSpace(N)
tspace = transform(space)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
F  = transformOp(space)

""" operators """
Â = diffusionOp(ν, tspace, discr)

v = 1.0;
vel = @. x*0 + v

Ĉ = advectionOp((F * vel,), tspace, discr)
F̂ = -Ĉ

Â = cache_operator(Â, im*k)
F̂ = cache_operator(F̂, im*k)

""" IC """
function uIC(x)
    @. sin(1x)
end
u0 = uIC(x)

""" time discr """
tspan = (0.0, 10.0)
#tsave = (0, π/4, π/2, 3π/4, 2π,)
tsave = 0:2:8 * pi
odealg = Tsit5()
prob = SplitODEProblem(Â, F̂, û0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave, abstol=1e-8, reltol=1e-8)

""" analysis """
pred = (F,) .\ sol.u
pred = hcat(pred...)

utrue(x,v,t) = uIC(@. x - v*t)
utr = utrue(x,v,sol.t[1])
for i=2:length(sol.u)
    ut = utrue(x,v,sol.t[i])
    global utr = hcat(utr, ut)
end

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, pred[:,i], legend=true)
end
display(plt)

err = norm(pred .- utr,Inf)
display(err)
@test err < 1e-8
#
