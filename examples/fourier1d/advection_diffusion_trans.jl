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
p = nothing

""" space discr """
space = FourierSpace(N)
tspace = transform(space)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
F  = transformOp(space)

""" operators """
v = 1.0;
vel = @. x*0 + v
vels = (F * vel,)

Â = diffusionOp(ν, tspace, discr)
Ĉ = advectionOp(vels, tspace, discr)
F̂ = NullOperator(tspace)
Dt = cache_operator(Â-Ĉ+F̂, im*k)

""" IC """
function uIC(x)
    @. sin(2x)
end
u0 = uIC(x)
û0 = F * u0

""" time discr """
tsave = 0: pi/4: 4pi
tspan = (tsave[begin], tsave[end])
odealg = Tsit5()
prob = ODEProblem(Dt, û0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave, abstol=1e-10, reltol=1e-10)

""" analysis """
pred = (F,) .\ sol.u
pred = hcat(pred...)

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, pred[:,i], legend=true)
end
display(plt)
#
