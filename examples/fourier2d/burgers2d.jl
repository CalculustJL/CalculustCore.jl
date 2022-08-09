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

using OrdinaryDiffEq, Plots
using ComponentArrays, CUDA

T = Float32
nx = ny = 128
N = nx * ny
ν = 1e-3 |> T
p = nothing

Nmodes = 16
tol = 1e-6 |> T

odealg = Tsit5()
odealg = SSPRK43()

""" spatial discr """
space = FourierSpace(nx, ny) |> T
discr = Collocation()
x, y = points(space)

odecb = begin
    function affect!(int)
        println(
                "[$(int.iter)] \t Time $(round(int.t; digits=8))s"
               )
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

""" IC """
u0 = begin
    X = truncationOp(space, (Nmodes//nx, Nmodes//ny))
    vx0 = X * rand(T, size(x)...)
    vy0 = X * rand(T, size(x)...)

    ComponentArray(;vx=vx0, vy=vy0)
end

ps = ComponentArray(vel=u0)
space = make_transform(space, u0.vx; p=ps)

# GPU
#CUDA.allowscalar(false)
#space = space |> gpu
#x, y = points(space)
#u0 = u0 |> gpu
#ps = ps |> gpu

""" spce ops """
Ax = diffusionOp(ν, space, discr)
Ay = diffusionOp(ν, space, discr)

Cx = advectionOp((zero(x), zero(x)), space, discr;
                 vel_update_funcs=(
                                   (v,u,p,t) -> copy!(v, p.vel.vx),
                                   (v,u,p,t) -> copy!(v, p.vel.vy),
                                  )
                )

Cy = advectionOp((zero(x), zero(x)), space, discr;
                 vel_update_funcs=(
                                   (v,u,p,t) -> copy!(v, p.vel.vx),
                                   (v,u,p,t) -> copy!(v, p.vel.vy),
                                  )
                )

Fx = NullOperator(space)
Fy = NullOperator(space)

Dtx = cache_operator(Ax-Cx+Fx, x)
Dty = cache_operator(Ay-Cy+Fy, x)

function ddt(du, u, p, t)
    ps = ComponentArray(vel=u)

    Dtx(du.vx, u.vx, ps, t)
    Dty(du.vy, u.vy, ps, t)

    du
end

""" time discr """
tspan = (0.0, 10.0) .|> T
tsave = range(tspan...; length=100)
prob  = ODEProblem(ddt, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave, abstol=1f-2, reltol=1f-2, callback=odecb)
@time sol = solve(prob, odealg, saveat=tsave, abstol=tol, reltol=tol, callback=odecb)
@time sol = solve(prob, odealg, saveat=tsave, abstol=tol, reltol=tol)
@time sol = solve(prob, odealg, saveat=tsave, abstol=tol, reltol=tol)

pred = Array(sol)
vx = @views pred[:vx, :]
vy = @views pred[:vy, :]

anim = animate(vx, space, sol.t)
filename = joinpath(dirname(@__FILE__), "burgers_x" * ".gif")
gif(anim, filename, fps=20)

anim = animate(vy, space, sol.t)
filename = joinpath(dirname(@__FILE__), "burgers_y" * ".gif")
gif(anim, filename, fps=20)
#
