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
using ComponentArrays, UnPack

nx = 32
ny = 32
ν = 1e-2
p = nothing

""" spatial discr """
space = FourierSpace(nx, ny)
discr = Collocation()

x, y = points(space)

Ax = diffusionOp(ν, space, discr)
Ay = diffusionOp(ν, space, discr)

Cx = advectionOp((zero(x), zero(x)), space, discr;
                 vel_update_funcs=(
                                   (v,u,p,t; vel=nothing) -> copy!(v, vel.vx),
                                   (v,u,p,t; vel=nothing) -> copy!(v, vel.vy),
                                  )
                )

Cy = advectionOp((zero(x), zero(x)), space, discr;
                 vel_update_funcs=(
                                   (v,u,p,t; vel=nothing) -> copy!(v, vel.vx),
                                   (v,u,p,t; vel=nothing) -> copy!(v, vel.vy),
                                  )
                )

Fx = -Cx
Fy = -Cy

Dtx = cache_operator(Ax+Fx, x)
Dty = cache_operator(Ay+Fy, x)

""" IC """
u0 = begin
    vx0 = @. sin(2x) * sin(2y)
    vy0 = @. sin(3x) * sin(3y)

    ComponentArray(vx=vx0, vy=vy0)
end

function ddt(du, u, p, t)
    SciMLOperators.update_coefficients!(Dtx, u.vx, p, t; vel=u)
    SciMLOperators.update_coefficients!(Dtx, u.vy, p, t; vel=u)

    mul!(du.vx, Dtx, u.vx)
    mul!(du.vy, Dty, u.vy)

    du
end

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=10)
odealg = Tsit5()
prob = ODEProblem(ddt, u0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

pred = Array(sol)
vx = @views pred[:vx, :]
vy = @views pred[:vx, :]

anim = animate(vx, space)
filename = joinpath(dirname(@__FILE__), "burgers_x" * ".gif")
gif(anim, filename, fps=5)

anim = animate(vy, space)
filename = joinpath(dirname(@__FILE__), "burgers_y" * ".gif")
gif(anim, filename, fps=5)
#
