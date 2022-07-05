#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, LinearAlgebra, Random
using DiffEqGPU
using Plots

using CUDA
CUDA.allowscalar(false)

N = 1024
ν = 1f-3
p = nothing

Random.seed!(0)
function uIC(space)
    x = points(space)[1]
    X = truncationOp(space,1//20)

    u0 = if x isa CUDA.CuArray
        X * CUDA.rand(size(x)...)
    else
        X * rand(size(x)...)
    end

    u0
end

function solve_burgers(N, ν, p;
                       uIC=uIC,
                       tspan=(0f0, 10f0),
                       nsims=1,
                       nsave=100,
                       odealg=SSPRK43(),
#                      odealg=GPUTsit5(),
                      )

    """ space discr """
    space = FourierSpace(N) |> gpu
    discr = Collocation()

    (x,) = points(space)

    """ IC """
    u0 = uIC(space)
    u0 = u0 * gpu(ones(1,nsims))
    space = make_transform(space, u0; p=p)

    """ operators """
    A = diffusionOp(ν, space, discr)

    function burgers!(v, u, p, t)
        copyto!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
    end

    C = advectionOp((zero(u0),), space, discr; vel_update_funcs=(burgers!,))
    F = -C + forcingOp(zero(u0), space, discr; f_update_func=forcing!)

    A = cache_operator(A, u0)
    F = cache_operator(F, u0)

    """ time discr """
#   odefunc = SplitFunction(A, F)
    odefunc = cache_operator(A+F, u0)

    tsave = range(tspan...; length=nsave)
    prob = ODEProblem(odefunc, u0, tspan, p; reltol=1f-6, abstol=1f-6)
    @time sol = solve(prob, odealg, saveat=tsave)

    sol, space
end

function plot_sol(pred, time, x)
    plt = plot()
    for i=1:length(time)
        plot!(plt, x, pred[:,i], legend=false)
    end
    plt
end

function anim8(pred, time, x)
    ylims = begin
        u = @views pred[:,1]
        mi = minimum(u)
        ma = maximum(u)
        buf = (ma-mi)/5
        (mi-buf, ma+buf)
    end
    anim = @animate for i=1:length(time)
        plt = plot(x, pred[:,i], legend=false, ylims=ylims)
    end
end

sol, space = solve_burgers(N, ν, p)

pred = Array(sol) |> cpu
time = sol.t |> cpu
(x,) = points(space) |> cpu

#plt = plot_sol(pred, time, x)
#display(plt)
#anim = anim8(pred, time, x)
#gif(anim, "examples/fourier/a.gif", fps= 20)
#display(plt)
nothing
#
