#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearAlgebra, ComponentArrays
using Lux, Random
using SciMLSensitivity, Zygote
using Optimization, OptimizationOptimJL, OptimizationOptimisers
using Plots

N = 128
ν = 1f-1
odealg = SSPRK43()

""" NN """
model = Lux.Chain(
                  Lux.Dense(1, 10),
                  Lux.Dense(10, 1),
                 )

rng = Random.default_rng()
ps, st = Lux.setup(rng, model)
ps = ComponentArray(ps)

""" space discr """
space  = FourierSpace(N) |> Float32
space  = make_transform(space; p=ps)
discr  = Collocation()

(x,) = points(space)

A = diffusionOp(ν, space, discr)
burgers!(v, u, p, t) = copy!(v, u)
forcing!(f, u, p, t) = lmul!(false, f)
C = advectionOp((zero(x),), space, discr; vel_update_funcs=(burgers!,))
F = -C + forcingOp(zero(x), space, discr; f_update_func=forcing!)

Dt = cache_operator(A+F, x)

u0 = @. sin(x)
tspan = (0f0, 10f0)
tsteps = range(tspan..., length=10)

function dudt(u, p, t; Dt=Dt, space=space, st=st)
    Zygote.ignore() do
        SciMLOperators.update_coefficients!(Dt, u, p, t)
    end

    du1 = Dt * u

    x = points(space)[1]
    du2 = model(x', p, st)[1] |> vec

    du1 + du2
end

prob = ODEProblem(dudt, u0, tspan, saveat=tsteps)
sense = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))

function predict(ps; prob=prob, odealg=odealg, sense=sense)
    solve(prob, odealg, p=ps, sensealg=sense) |> Array
end

function loss(p)
    pred = predict(p)
    loss = sum(abs.(pred .- 1.0))

    loss, pred
end

function cb(p, l, pred; doplot=false, space=space)
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

# dummy
println("fwd"); cb(ps,loss(ps)...;doplot=true)
println("bwd"); Zygote.gradient(p -> loss(p)[1], ps) |> display

""" optimization """
adtype = Optimization.AutoZygote()
# x=object to optimize
# p=parameters for optimization loop
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ps)

optres = Optimization.solve(
                            optprob,
                            ADAM(0.05),
                            callback=cb,
                            maxiters=50,
                           )

optprob = remake(optprob,u0 = optres.u)

println("BFGS")
optres = Optimization.solve(optprob,
                            Optim.BFGS(initial_stepnorm=0.01),
                            callback=cb,
                            allow_f_increases = false,
                           )

cb(optres.u, loss(optres.u)...; doplot=true)
#
