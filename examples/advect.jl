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
#   @. sin(2x)
    u0 = rand(ComplexF64, size(k))
    u0[20:end] .= 0
    ftr \ u0
end
u0 = uIC(x)

""" time discr """
tspan = (0.0, 10.0)
tsave = (0, 2π)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(A, F, u0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """

@show norm(sol.u[2] .- sol.u[1], Inf)

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=true)
end
plt
#
