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

using CUDA
CUDA.allowscalar(false)

N = 128
ν = 1f-2
p = ()

""" space discr """
space = FourierSpace(N)
space = gpu(space)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
ftr  = transformOp(space)

α = 5
u0 = @. sin(α*x)

A = diffusionOp(ν, space, discr)
F = SciMLOperators.NullOperator(space)

A = cache_operator(A, x)
F = cache_operator(F, x)

function implicit(u, p, t;op=A)
    op * u
end

function explicit(u, p, t; op=F)
    op * u
end

""" time discr """
tspan = (0f0, 10f0) 
tsave = range(tspan...; length=10)
odealg = Tsit5()
#prob = SplitODEProblem(A, F, u0, tspan, p)
prob = SplitODEProblem(implicit, explicit, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave, reltol=1f-6)

""" analysis """
pred = cpu(Array(sol))

u0 = Array(u0)
utrue(t) = @. u0 * (exp(-ν*α^2*t))
ut = utrue(sol.t[1])
for i=2:length(sol.t)
    utt = utrue(sol.t[i])
    global ut = hcat(ut, utt)
end

plt = plot()
x = Array(x)
for i=1:length(sol.u)
    plot!(plt, x, pred[:,i], legend=false)
end
display(plt)

err = norm(pred .- ut, Inf)
display(err)
@test err < 1e-5
#
