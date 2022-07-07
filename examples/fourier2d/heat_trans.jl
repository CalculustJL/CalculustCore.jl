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
space  = FourierSpace(N)
tspace = transform(space)
discr  = Collocation()

(x,) = points(space)
(k,) = points(tspace)
iftr = transformOp(tspace)
ftr  = transformOp(space)

α = 5
u0 = @. sin(α*x)
û0 = ftr * u0

Â = diffusionOp(ν, tspace, discr)
F̂ = SciMLOperators.NullOperator(tspace)

Â = cache_operator(Â, k)
F̂ = cache_operator(F̂, k)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(Â, F̂, û0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
pred = iftr.(sol.u, nothing, 0)
pred = hcat(pred...)

utrue(t) = @. u0 * (exp(-ν*α^2*t))
ut = utrue(sol.t[1])
for i=2:length(sol.t)
    utt = utrue(sol.t[i])
    global ut = hcat(ut, utt)
end

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, pred[:,i], legend=false)
end
display(plt)

err = norm(pred .- ut, Inf)
@test err < 1e-6
#
