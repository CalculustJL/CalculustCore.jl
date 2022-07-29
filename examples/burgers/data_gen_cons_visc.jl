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
using CUDA, Random, JLD2, Plots

Random.seed!(0)
CUDA.allowscalar(false)

function uIC(space; truncation_frac=N_target/N)
    x = points(space)[1]
    X = truncationOp(space, (truncation_frac,))

    u0 = if x isa CUDA.CuArray
        X * CUDA.rand(size(x)...)
    else
        X * rand(size(x)...)
    end

    u0
end

odecb = begin
    function affect!(int)
        println(
                "[$(int.iter)] \t Time $(round(int.t; digits=8))s"
               )
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

function solve_burgers1D(N, ν, p;
                         uIC=uIC,
                         tspan=(0f0, 10f0),
                         nsims=10,
                         nsave=100,
                         odealg=SSPRK43(),
                        )

    """ space discr """
    space = FourierSpace(N) |> gpu
    discr = Collocation()

    (x,) = points(space)

    """ IC """
    u0 = [uIC(space) for i=1:nsims]
    u0 = hcat(u0...)
    space = make_transform(space, u0; p=p)

    """ operators """
    function burgers!(v, u, p, t)
        copyto!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
    end

    A = diffusionOp(ν, space, discr)
    C = advectionOp((zero(u0),), space, discr; vel_update_funcs=(burgers!,))
    F = forcingOp(zero(u0), space, discr; f_update_func=forcing!)
    odefunc = cache_operator(A-C+F, u0) |> ODEFunction

    """ time discr """
    tsave = range(tspan...; length=nsave)
    prob = ODEProblem(odefunc, u0, tspan, p; reltol=1f-6, abstol=1f-6)
    @time sol = solve(prob, odealg, saveat=tsave, callback=odecb)

    sol, space
end

function process(datafile; fps=20)
    data = jldopen(datafile)

    sp_coarse = data["sp_coarse"]
    sp_dense  = data["sp_dense"]
    u_coarse  = data["u_coarse"]
    u_dense   = data["u_dense"]
    t = data["t"]

    datadir   = dirname(datafile)
    densedir  = mkpath(joinpath(datadir, "dense"))
    coarsedir = mkpath(joinpath(datadir, "coarse"))

    for i=1:size(u_coarse, 2)
        namec = joinpath(coarsedir, "trajectory" * "$i")
        uc = @view u_coarse[:,i,:]
        animc = animate(uc, sp_coarse, t)
        gif(animc, namec * ".gif"; fps=fps)

        named = joinpath(densedir, "trajectory" * "$i")
        ud = @view u_dense[:,i,:]
        animd = animate(ud, sp_dense, t)
        gif(animd, named * ".gif"; fps=fps)
    end

    nothing
end

function datagen_burgers1D(N, ν, p, N_target, filename; kwargs...)

    sol, space = solve_burgers1D(N, ν, p; kwargs...)

    u_dense  = sol |> Array
    sp_dense = make_transform(space, u_dense) |> cpu

    sp_coarse = begin
        sz = (N_target, size(u_dense)[2:end]...)
        u  = similar(u_dense, sz)

        make_transform(FourierSpace(N_target), u)
    end

    J = interpOp(sp_coarse, sp_dense)

    u_coarse = J * u_dense
    t = sol.t |> cpu

    datadir = joinpath(@__DIR__, filename) |> mkpath
    datafile = joinpath(datadir, filename * ".jld2")

    jldsave(datafile; sp_coarse, sp_dense, u_coarse, u_dense, t)

    datafile
end

#########################
N = 1024
ν = 1f-3
p = nothing

N_target = 128

filename = "burgers_nu1em3_n1024"
datafile = datagen_burgers1D(N, ν, p, N_target, filename)
process(datafile)
#
