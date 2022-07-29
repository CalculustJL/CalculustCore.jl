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
1D Burgers + Closure equation

∂t(vx) = -vx * ∂x(vx) + ν∂xx(vx) + ∇⋅η
∂t(η) = -α*u∂x(η) + β*ν∂xx(η) + NN1(vx)

u(t=0) = u0 (from data)
η(t=0) = NN0(vx0)

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
    vx_data, t_data = ut_from_data(filename)
    vx0 = vx_data[:,:,1]
    n_data = length(vx0)

    """ initial conditions """
    u0 = ComponentArray(;
                        vx=vx0,
                        η=zero(vx0),
                       ) |> gpu

    """ operators """
    space = make_transform(space, u0.vx; isinplace=false, p=u0)

    Dx = gradientOp(space, discr)[1]
    Dx = cache_operator(Dx, u0.η)

    Ddt_vx = begin
        A = diffusionOp(ν, space, discr)

        function burgers!(v, u, p, t)
            copyto!(v, p.vx)
        end

        function forcing!(f, u, p, t)
            mul!(f, Dx, p.η)
        end

        C = advectionOp((zero(u0.vx),), space, discr; vel_update_funcs=(burgers!,))
        F = forcingOp(zero(u0.vx), space, discr; f_update_func=forcing!)

        cache_operator(A-C+F, u0.vx)
    end

    Ddt_η = begin
        A = diffusionOp(ν, space, discr)

        function vel!(v, u, p, t)
            copy!(v, p.u.vx)
        end

        C = advectionOp((zero(u0.η),), space, discr; vel_update_funcs=(vel!,))

        cache_operator(A-C, u0.η)
    end

    """ time discr """
    function dudt(u, p, t)
#       Zygote.ignore() do
#           SciMLBase.update_coefficients!(Dt_vx, u.vx, u, t)
#           SciMLBase.update_coefficients!(Dt_η , u.η , u, t)
#       end

#       dvx = Ddt_vx * u.vx
#       dη  = Ddt_η  * u.η

#       ComponentArray(;
#                      vx=dvx,
#                      η=dη,
#                     )
        zero(u)
    end

    tspan = (t_data[1], t_data[end])
    prob  = ODEProblem(dudt, u0, tspan, p; reltol=1f-6, abstol=1f-6)
    sense = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))

    function predict(p; callback=ode_cb)
#       prob = remake(prob,
#                     u0=ComponentArray(;
#                                       vx=,
#                                       η=,
#                                      )
#                    )
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

ode_cb = begin
    function affect!(int)
        println(int.t)
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

##############################################
name = "burgers_nu1em3_n1024"
filename = joinpath(@__DIR__, name * ".jld2")

N = 128
ν = 1f-3

model, ps, st = begin
    nn_η_init = Lux.Chain(
                          Lux.Dense(N, N),
                          Lux.Dense(N, N),
                         )

    p_η_init, st_η_init = Lux.setup(rng, nn_η_init)

    nn_η_forcing = Lux.Chain(
                             Lux.Dense(N,N),
                             Lux.Dense(N,N),
                            )

    p_η_forcing, st_η_forcing = Lux.setup(rng, nn_η_forcing)

    α = rand(Float32)
    β = rand(Float32)

    model = (
             η_init = nn_η_init,
             η_forcing = nn_η_forcing,
            )

    ps = ComponentArray(
                        η_init=p_η_init,
                        η_forcing=p_η_forcing,
                        α=α,
                        β=β,
                       ) |> gpu

    st = (;
          η_init = st_η_init,
          η_forcing = st_η_forcing,
         )

    model, ps, st
end

predict, loss, space = setup_burgers1d(N, ν, filename; p=ps, model=model);

# dummy calls
println("fwd"); @time opt_cb(ps, loss(ps)...;doplot=false)
#println("bwd"); @time Zygote.gradient(p -> loss(p)[1], ps) |> display

#ps = train(loss, ps)
#