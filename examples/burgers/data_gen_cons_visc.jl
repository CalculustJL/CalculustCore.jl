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
using CUDA, Random, JLD2

Random.seed!(0)
CUDA.allowscalar(false)

N = 8192
ν = 1f-5
p = nothing

N = 1024
ν = 1f-3
p = nothing

N_target = 128

name = "burgers_nu1em5_n8192"

function uIC(space; truncation_frac=N_target/N)
    x = points(space)[1]
    X = truncationOp(space, truncation_frac)

    u0 = if x isa CUDA.CuArray
        X * CUDA.rand(size(x)...)
    else
        X * rand(size(x)...)
    end

    u0
end

function solve_burgers(N, ν, p;
                       uIC=uIC,
                       tspan=(0f0, 100f0),
                       nsims=100,
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

sol, space = solve_burgers(N, ν, p)

pred = Array(sol) |> cpu # [npts, nsims, nsave]
time = sol.t |> cpu
(x,) = points(space) |> cpu

filename = joinpath(dirname(@__FILE__), name * ".jld2")
jldsave(filename; pred, time, space)

nothing
#
