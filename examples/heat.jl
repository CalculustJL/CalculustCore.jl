#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces
using OrdinaryDiffEq, LinearSolve, LinearAlgebra
using Plots

N = 1024
ν = 1e-2
p = ()

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)
tr = space.transforms
k = modes(space)

α = 5
u0 = @. sin(α*x)

A = diffusionOp(ν, space, discr)
F = SciMLOperators.NullOperator(space)

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
pred = Array(sol)

utrue(t) = @. u0 * (exp(-ν*α^2*t))
ut = utrue(sol.t[1])
for i=2:length(sol.t)
    utt = utrue(sol.t[i])
    global ut = hcat(ut, utt)
end

norm(pred .- ut, Inf) |> display

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=false)
end
display(plt)

nothing
