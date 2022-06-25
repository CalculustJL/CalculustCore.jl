#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, LinearAlgebra, Sundials, Random
using Plots

N = 1024
ν = 1e-3
p = ()

Random.seed!(0)
function uIC(x, ftr, k)
#   u0 = @. sin(2x) + sin(3x) + sin(5x)
#   u0 = @. sin(x-π)

    u0 = begin
        u  = 2*rand(size(x)...)
        uh = ftr * u
        uh[20:end] .= 0
        ftr \ uh
    end

    u0
end

function solve_burgers(N, ν, p;
                       uIC=uIC,
                       tspan=(0.0, 10.0),
                       nsave=100,
                      )

    """ space discr """
    space = FourierSpace(N)
    discr = Collocation()

    (x,) = points(space)
    (k,) = modes(space)
    ftr  = transforms(space)

    """ IC """
    u0 = uIC(x, ftr, k)

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
    odealg = CVODE_BDF(method=:Functional)
    tsave = range(tspan...; length=nsave)
    prob = SplitODEProblem(A, F, u0, tspan, p)
    @time sol = solve(prob, odealg, saveat=tsave)

    sol, space
end

function plot_sol(sol::ODESolution, space::FourierSpace)
    x = points(space)[1]
    plt = plot()
    for i=1:length(sol)
        plot!(plt, x, sol.u[i], legend=false)
    end
    plt
end

function anim8(sol::ODESolution, space::FourierSpace)
    x = points(space)[1]
    ylims = begin
        u = sol.u[1]
        mi = minimum(u)
        ma = maximum(u)
        buf = (ma-mi)/3
        (mi-buf, ma+buf)
    end
    anim = @animate for i=1:length(sol)
        plt = plot(x, sol.u[i], legend=false, ylims=ylims)
    end
end

sol, space = solve_burgers(N, ν, p)
#plt = plot_sol(sol, space)
anim = anim8(sol, space)
gif(anim, "examples/fourier/a.gif", fps= 20)

#
nothing
