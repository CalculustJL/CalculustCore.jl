#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearAlgebra
using Test, Plots

nx = 32
ny = 32
ν = 0e0
p = nothing

""" space discr """
space = FourierSpace(nx, ny)
discr = Collocation()

(x,y) = points(space)

""" operators """
A = diffusionOp(ν, space, discr)

vx = 1.0; velx = @. x*0 + vx
vy = 1.0; vely = @. x*0 + vy
C = advectionOp((velx, vely), space, discr)
F = -C

A = cache_operator(A, x)
F = cache_operator(F, x)

""" IC """
uIC(x,y) = @. sin(1x) * sin(1y)
u0 = uIC(x,y)

""" time discr """
tspan = (0.0, 10.0)
tsave = (0, π/4, π/2, 3π/4, 2π,)
odealg = Tsit5()
prob = SplitODEProblem(A, F, u0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave, abstol=1e-8, reltol=1e-8)

""" analysis """
pred = Array(sol)

utrue(x, y, vx, vy, t) = uIC(x .- vx*t, y .- vy*t)
utr = utrue(x,y,vx,vy,sol.t[1])
for i=2:length(sol.u)
    ut = utrue(x, y, vx, vy, sol.t[i])
    global utr = hcat(utr, ut)
end

anim = animate(pred, space)
filename = joinpath(dirname(@__FILE__), "advect" * ".gif")
gif(anim, filename, fps=5)

err = norm(pred .- utr,Inf)
@test err < 1e-8
#
