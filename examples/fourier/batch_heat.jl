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
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)

α = 5
u0 = @. sin(α*x)

u0 = u0 * ones(1,2)
space = make_transform(space, u0)

A = diffusionOp(ν, space, discr)
F = SciMLOperators.NullOperator(space)

A = cache_operator(A, u0)
F = cache_operator(F, u0)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
pred = Array(sol)

ut = similar(pred)
utrue(t) = @. u0 * (exp(-ν*α^2*t))

for i=1:length(sol)
    ut[:,:,i] = utrue(sol.t[i])
end

err = norm(pred .- ut, Inf)
@test err < 1e-6
#
