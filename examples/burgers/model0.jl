#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, CUDA, LinearAlgebra, ComponentArrays
using Lux, Random, JLD2, SciMLSensitivity, Zygote
using Optimization, OptimizationOptimJL, OptimizationOptimisers, Optimisers
using Plots

CUDA.allowscalar(false)

rng = Random.default_rng()
Random.seed!(rng, 0)

"""
1D Burgers + NN forcing

∂t(u) + u∂x(u) = νΔ(u) + NN(u)
"""

""" data """
function ut_from_data(filename)
    data = jldopen(filename)
    
    t = data["t"]
    u = data["u_coarse"]

    u, t
end

""" space discr """
function setup_burgers1d(N, ν, filename;
                         p=nothing,
                         model=nothing,
                         odealg=SSPRK43(),
                         ode_cb=nothing,
                        )

    model = model isa Nothing ? (u, p, t, space) -> zero(u) : model

    """ space discr """
    space = FourierSpace(N) |> gpu
    discr = Collocation()

    (x,) = points(space)

    """ get data """
    u_data, t_data = ut_from_data(filename)
    u_data = gpu(u_data)
    u0 = @views u_data[:,:,1]
    n_data = length(u_data)

    """ operators """
    space = make_transform(space, u0; isinplace=false, p=p)
    A = diffusionOp(ν, space, discr)

    function burgers!(v, u, p, t)
        copyto!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
    end

    C = advectionOp((zero(u0),), space, discr; vel_update_funcs=(burgers!,))
    F = -C + forcingOp(zero(u0), space, discr; f_update_func=forcing!)

    Dt = cache_operator(A+F, u0)

    """ time discr """
    function dudt(u, p, t)
        Zygote.ignore() do
            SciMLBase.update_coefficients!(Dt, u, p, t)
        end

        du1 = Dt * u
        du2 = model(u, p, t, space)

        du1 + du2
    end

    tspan = (t_data[1], t_data[end])
    prob = ODEProblem(dudt, u0, tspan, p; reltol=1f-6, abstol=1f-6)
    sense = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))

    function predict(p; callback=ode_cb)
        solve(prob,
              odealg,
              p=p,
              sensealg=sense,
              callback=callback,
              saveat=t_data,
             ) |> CuArray
    end

    function loss(p)
        pred = predict(p)
        loss = sum(abs2.(u_data .- pred)) / n_data

        loss, pred
    end

    predict, loss, space
end

function train(loss, p;
               optalg=Optimisers.ADAM(1f-3),
               niters=100,
              )

    opt_f  = p -> loss(p)[1]
    opt_st = Optimisers.setup(optalg, p)

    stime = time()
    for iter in 1:niters

        # loss, gradient
        l, back = pullback(opt_f, p)
        g = back(one(l))[1]

        # update 
        opt_st, p = Optimisers.update(opt_st, p, g)

        # logging
        ttime = time() - stime
        opt_cb(p, l, pred; steptime=ttime)

    end

    p
end

ode_cb = begin
    function affect!(int)
        println(int.t)
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

function opt_cb(p, l, pred; doplot=false, space=space, steptime=nothing)
    println(
            "[$iter/$niters] \t Time $(round(ttime; digits=2))s \t Loss: " *
            "$(round(l; digits=8)) \t "
           )

    return false
end

##############################################
name = "burgers_nu1em3_n1024"
filename = joinpath(@__DIR__, name * ".jld2")

N = 128
ν = 1f-3

""" NN """
model, ps, st = begin
    w = N
    nn = Lux.Chain(
                   Lux.Dense(w, w),
                   Lux.Dense(w, w),
                  )

    rng = Random.default_rng()
    Random.seed!(rng, 0)

    ps, st = Lux.setup(rng, nn)
    ps = ComponentArray(ps) |> gpu
    st = st |> gpu

    function model(u, p, t, space)
        nn(u, p, st)[1]
    end

    model, ps, st
end

predict, loss, space = setup_burgers1d(N, ν, filename; p=ps, model=model);

# dummy calls
println("fwd"); @time opt_cb(ps, loss(ps)...;doplot=false)
println("bwd"); @time Zygote.gradient(p -> loss(p)[1], ps) |> display

#ps = train(loss, ps)
#
