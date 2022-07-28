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
function ut_from_data(datafile)
    data = jldopen(datafile)
    
    t = data["t"]
    u = data["u_coarse"]

    u, t
end

odecb = begin
    function affect!(int)
        println(int.t)
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

function optcb(p, l, pred;
                doplot=false,
                space=space,
                steptime=nothing,
                iter=nothing,
                niter=nothing,
               )

    steptime = steptime isa Nothing ? 0.0 : steptime
    iter = iter isa Nothing ? 0 : iter
    niter = niter isa Nothing ? 0 : niter

    println(
            "[$iter/$niter] \t Time $(round(steptime; digits=2))s \t Loss: " *
            "$(round(l; digits=8)) \t "
           )

    return false
end

""" space discr """
function setup_burgers1d(N, ν, datafile;
                         p=nothing,
                         model=nothing,
                         odealg=SSPRK43(),
                         odecb=nothing,
                        )

    model = model isa Nothing ? (u, p, t, space) -> zero(u) : model

    """ space discr """
    space = FourierSpace(N) |> gpu
    discr = Collocation()

    (x,) = points(space)

    """ get data """
    u_data, t_data = ut_from_data(datafile)
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
        du2 = 1f-4*model(u, p, t, space)

        du1 + du2
    end

    tspan = (t_data[1], t_data[end])
    prob = ODEProblem(dudt, u0, tspan, p; reltol=1f-4, abstol=1f-4)
    sense = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))

    function predict(p; callback=odecb)
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
               alg=Optimisers.Adam(1f-2),
               maxiters=1000,
               callback=optcb,
              )

    adtype = Optimization.AutoZygote()
    # x=object to optimize
    # p=parameters for optimization loop
    optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, ps)

    optres = Optimization.solve(optprob, alg, callback=callback, maxiters=maxiters)

    optres.u
end

##############################################
name = "burgers_nu1em3_n1024"
datafile = joinpath(@__DIR__, name * ".jld2")
savefile = joinpath(@__DIR__, "model0" * ".jld2")

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

#odealg = Tsit5()
odealg = SSPRK43()
predict, loss, space = setup_burgers1d(N, ν, datafile; odealg=odealg, p=ps, model=model);

# dummy calls
println("fwd"); @time optcb(ps, loss(0*ps)...;doplot=false)
println("bwd"); @time Zygote.gradient(p -> loss(p)[1], 0*ps) |> display

ps = train(loss, ps; alg=ADAM(1f-1), maxiters=100)
ps = train(loss, ps; alg=ADAM(1f-2), maxiters=1000)
#ps = train(loss, ps; alg=ADAM(1f-3), maxiters=5000)

l, pred = loss(ps)
display(l)

pred = cpu(pred)
space = cpu(space)

jldsave(savefile; ps)

# plot
for i=1:10
    name = "trajectory" * "$i"
    u = @view pred[:,i,:]
    anim = animate(u, space)
    gif(anim, name * ".gif"; fps=10)
end

#
