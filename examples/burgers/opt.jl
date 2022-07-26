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
using Optimization, OptimizationOptimJL, OptimizationOptimisers
using Plots

Random.seed!(0)
CUDA.allowscalar(false)

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

function opt_cb(p, l, pred; doplot=false, space=space)
    println(l)

    if doplot
        plt = plot()
        for i=1:size(pred,2)
            x = points(space)[1]
            plot!(plt, x, pred[:,i])
        end
        display(plt)
    end
    return false
end

##############################################
name = "burgers_nu1em3_n1024"
filename = joinpath(@__DIR__, name * ".jld2")

N = 128
ν = 1f-3

""" NN """
model, ps = begin
    w = N
    nn = Lux.Chain(
                   Lux.Dense(w, w),
                   Lux.Dense(w, w),
                  )

    rng = Random.default_rng()
    ps, st = Lux.setup(rng, nn)
    ps = ComponentArray(ps) |> gpu

    function model(u, p, t, space)

        nn(u, p, st)[1]
    end

    model, ps
end

ode_cb = begin
    function affect!(int)
        println(int.t)
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

predict, loss, space = setup_burgers1d(N, ν, filename; p=ps, model=model);

## dummy
println("fwd"); @time opt_cb(ps, loss(ps)...;doplot=false)
println("bwd"); @time Zygote.gradient(p -> loss(p)[1], ps) |> display

#""" optimization """
#adtype = Optimization.AutoZygote()
## x=object to optimize
## p=parameters for optimization loop
#optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
#optprob = Optimization.OptimizationProblem(optf, ps)
#
#optres = Optimization.solve(
#                            optprob,
#                            ADAM(0.05),
#                            callback=cb,
#                            maxiters=50,
#                           )
#
#optprob = remake(optprob,u0 = optres.u)
#
#println("BFGS")
#optres = Optimization.solve(optprob,
#                            Optim.BFGS(initial_stepnorm=0.01),
#                            callback=cb,
#                            allow_f_increases = false,
#                           )
#
#cb(optres.u, loss(optres.u)...; doplot=true)
#
