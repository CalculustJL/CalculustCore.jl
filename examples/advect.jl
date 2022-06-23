#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces, LinearAlgebra
using OrdinaryDiffEq, LinearSolve
using Plots

N = 1024
ν = 0e0
p = ()

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)
ftr = space.transforms
k = modes(space)

""" operators """
A = diffusionOp(ν, space, discr)

v = 1.0; vel     = @. x*0 + v
f = 0.0; forcing = @. x*0 + f
C = advectionOp((vel,), space, discr)
F = AffineOperator(C, forcing)

A = cache_operator(A, x)
F = cache_operator(F, x)

""" IC """
function uIC(x)
    @. sin(2x)
end

u0 = uIC(x)
function uT(t, v)
    xx = @. x - v*t
    uIC(xx)
end

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
utrue = uT(sol.t[1], v)
for i=2:length(sol.t)
    ut = uT(sol.t[i], v)
    global utrue = hcat(utrue, ut)
end

pred = Array(sol)
@show norm(pred - utrue)

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=false)
end
plt
