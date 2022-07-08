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

nx = 32
ny = 32
ν = 1e-0
p = nothing

""" space discr """
space = FourierSpace(nx, ny)
discr = Collocation()

x, y = points(space)

u0 = @. sin(2x) * sin(2y)

A = diffusionOp(ν, space, discr)

velx = @. x*0 + 1.0
vely = @. x*0 + 1.0
C = advectionOp((velx, vely), space, discr)
F = -C

A = cache_operator(A, x)
F = cache_operator(F, x)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Tsit5()
prob = SplitODEProblem(A, F, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

anim = animate(pred, space)
filename = joinpath(dirname(@__FILE__), "advection_diffusion" * ".gif")
gif(anim, filename, fps=5)

