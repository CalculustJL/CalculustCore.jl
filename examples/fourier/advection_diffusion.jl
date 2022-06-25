#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve
using Plots

N = 1024
ν = 1e-2
p = ()

""" space discr """
space = FourierSpace(N)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
ftr  = transforms(space)

u0 = @. sin(2x)

#u0 = rand(ComplexF64, size(k))
#u0[20:end] .= 0
#u0 = ftr \ u0

A = diffusionOp(ν, space, discr)

v = @. x*0 + 1
C = advectionOp((v,), space, discr)
F = -C

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=false)
end
plt
