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
using Test, Plots

N = 128
ν = 0e0
p = ()

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
ftr = transformOp(space)

""" operators """
A = diffusionOp(ν, space, discr)

v = 1.0; vel = @. x*0 + v
C = advectionOp((vel,), space, discr)
F = -C

A = cache_operator(A, x)
F = cache_operator(F, x)

""" IC """
function uIC(x)
    @. sin(1x)
#   u0 = rand(ComplexF64, size(k))
#   u0[20:end] .= 0
#   ftr \ u0
end
u0 = uIC(x)

""" time discr """
tspan = (0.0, 10.0)
tsave = (0, π/4, π/2, 3π/4, 2π,)
#tsave = (0, 2π,)
odealg = Rodas5(autodiff=false)
#odealg = Tsit5()
prob = SplitODEProblem(A, F, u0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave, abstol=1e-8, reltol=1e-8)

""" analysis """
pred = Array(sol)

utrue(x,v,t) = uIC(@. x - v*t)
utr = utrue(x,v,sol.t[1])
for i=2:length(sol.u)
    ut = utrue(x,v,sol.t[i])
    global utr = hcat(utr, ut)
end

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=true)
    plot!(plt, x, utr[:,i], legend=true)
end
display(plt)

err = norm(pred .- utr,Inf)
@test err < 1e-8
#
