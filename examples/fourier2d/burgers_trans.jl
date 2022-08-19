#
"""
Solve the 2D Burgers equation

∂t(vx) = -(vx*∂x(vx) + vy*∂y(vx)) + ν*Δvx
∂t(vy) = -(vx*∂x(vy) + vy*∂y(vy)) + ν*Δvy
"""

using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, LinearAlgebra
using ComponentArrays

T = Float32
nx = ny = 128
N = nx * ny
ν = 1e-3 |> T
p = nothing

Nmodes = 16
tol = 1e-6 |> T

odealg = Tsit5()
odealg = SSPRK43()

""" space discr """
space  = FourierSpace(nx, ny)
discr  = Collocation()

x, y = points(space)
F = transformOp(space)

odecb = begin
    function affect!(int)
        println(
                "[$(int.iter)] \t Time $(round(int.t; digits=8))s"
               )
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

u0 = begin
    X = truncationOp(space, (Nmodes//nx, Nmodes//ny))
    vx0 = X * rand(T, size(x)...)
    vy0 = X * rand(T, size(x)...)

    ComponentArray(;vx=vx0, vy=vy0)
end

û0 = begin
    v̂x0 = F * u0.vx
    v̂y0 = F * u0.vy

    ComponentArray(;v̂x=v̂x0, v̂y=v̂y0)
end

ps = ComponentArray(vel=û0)
space = make_transform(space, F \ û0.v̂x; p=ps)
tspace = transform(space)

""" spce ops """
D̂tx, D̂ty = begin
    ca = û0.v̂x

    Âx = diffusionOp(ν, tspace, discr)
    Ây = diffusionOp(ν, tspace, discr)

    Ĉx = advectionOp((zero(ca), zero(ca)), tspace, discr;
                     vel_update_funcs=(
                                       (v̂,û,p,t) -> copy!(v̂, p.vel.v̂x),
                                       (v̂,û,p,t) -> copy!(v̂, p.vel.v̂y),
                                      )
                    )

    Ĉy = advectionOp((zero(ca), zero(ca)), tspace, discr;
                     vel_update_funcs=(
                                       (v̂,û,p,t) -> copy!(v̂, p.vel.v̂x),
                                       (v̂,û,p,t) -> copy!(v̂, p.vel.v̂y),
                                      )
                    )

    F̂x = NullOperator(tspace)
    F̂y = NullOperator(tspace)

    D̂tx = cache_operator(Âx-Ĉx+F̂x, ca)
    D̂ty = cache_operator(Ây-Ĉy+F̂y, ca)

    D̂tx, D̂ty
end

function ddt(dû, û, p, t)
    ps = ComponentArray(vel=û)

    D̂tx(dû.v̂x, û.v̂x, ps, t)
    D̂ty(dû.v̂y, û.v̂y, ps, t)

    dû
end

""" time discr """
tspan = (0.0, 10.0) .|> T
tsave = range(tspan...; length=100)
prob  = ODEProblem(ddt, û0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave, abstol=1f-2, reltol=1f-2, callback=odecb)
@time sol = solve(prob, odealg, saveat=tsave, abstol=tol, reltol=tol, callback=odecb)
@time sol = solve(prob, odealg, saveat=tsave, abstol=tol, reltol=tol)
@time sol = solve(prob, odealg, saveat=tsave, abstol=tol, reltol=tol)

pred = Array(sol)
v̂x = @views pred[:v̂x, :]
v̂y = @views pred[:v̂y, :]

space = make_transform(space, zeros(N, size(pred, 2)))
F = transformOp(space)
vx = F \ v̂x
vy = F \ v̂y

anim = animate(vx, space, sol.t)
filename = joinpath(dirname(@__FILE__), "burgers_trans_x" * ".gif")
gif(anim, filename, fps=20)

anim = animate(vy, space, sol.t)
filename = joinpath(dirname(@__FILE__), "burgers_trans_y" * ".gif")
gif(anim, filename, fps=20)
#
