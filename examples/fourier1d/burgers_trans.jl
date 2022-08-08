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

N = 1024
ν = 1e-2
p = nothing

""" space discr """
space  = FourierSpace(N)
tspace = transform(space)
discr  = Collocation()

(x,) = points(space)
(k,) = points(tspace)
F    = transformOp(space)

""" initial condition """
function uIC(space)
    x = points(space)[1]
    X = truncationOp(space, (32/N,))

    u0 = X * rand(size(x)...)
end
u0 = uIC(space)
û0 = F * u0

function burgers!(v, u, p, t)
    copy!(v, u)
end

function forcing!(v, u, p, t)
    lmul!(false, v)
end

Â = diffusionOp(ν, tspace, discr)
Ĉ = advectionOp((zero(û0),), tspace, discr; vel_update_funcs=(burgers!,))
F̂ = forcingOp(zero(û0), tspace, discr; f_update_func=forcing!)

Dt = cache_operator(Â-Ĉ+F̂, û0)

""" time discr """
tspan = (0.0, 10.0)
tsave = range(tspan...; length=100)
odealg = Rodas5(autodiff=false)
prob = ODEProblem(Dt, û0, tspan, p)

@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
pred = [F,] .\ sol.u
pred = hcat(pred...)

anim = animate(pred, space)
gif(anim, "burg_trans.gif", fps=20)
#
