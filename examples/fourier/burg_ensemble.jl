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
using Plots

N = 1024
ν = 1f-3
p = ()

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
                       nsave=100,
                       odealg=SSPRK43(),
                      )

    """ space discr """
    space = FourierSpace(N) |> gpu
    discr = Collocation()

    (x,) = points(space)

    """ IC """
    u0 = [uIC(space) for i=1:2]

    """ operators """
    A = diffusionOp(ν, space, discr)

    function burgers!(v, u, p, t)
        copy!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
#       f .= 1e-2*rand(length(f))
    end

    C = advectionOp((zero(x),), space, discr; vel_update_funcs=(burgers!,))
    F = -C + forcingOp(zero(x), space, discr; f_update_func=forcing!)

    A = cache_operator(A, x)
    F = cache_operator(F, x)

    """ time discr """
    odefunc = SplitFunction(A, F)

    tsave = range(tspan...; length=nsave)
    odeprob = ODEProblem(odefunc, u0[1], tspan, p; reltol=1f-8, abstol=1f-8)
    @time sol = solve(odeprob, odealg, saveat=tsave)

    # problems selector function
    function prob_func(odeprob, i, repeat)
        odeprob = remake(ode_prob, u0=u0[i])
    end
    
    eprob = EnsembleProblem(odeprob, prob_func = prob_func)
    esol  = solve(eprob,
                  odealg,
                  EnsembleGPUArray();
                  trajectories=length(u0),
                  saveat=tsave,
                 )

    esol, space
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

plt = plot_sol(pred, time, x)
display(plt)
anim = anim8(pred, time, x)
gif(anim, "examples/fourier/a.gif", fps= 20)
display(plt)
#
