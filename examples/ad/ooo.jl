#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearAlgebra, CUDA
using Lux, Random, ComponentArrays
using SciMLSensitivity, Zygote
using Optimization, OptimizationOptimJL, OptimizationOptimisers
using Plots

Random.seed!(0)
CUDA.allowscalar(false)

N = 128
odealg = SSPRK43()

model = Lux.Chain(
                  Lux.Dense(N, N),
                  Lux.Dense(N, N),
                 )

rng = Random.default_rng()
ps, st = Lux.setup(rng, model)
ps = ComponentArray(ps) |> gpu
st = st |> gpu

""" space discr """
space  = FourierSpace(N)
(x,) = points(space)
u0 = sin.(10x) |> gpu
tspan = (0f0, 1f0)
tsteps = range(tspan..., length=10)

function dudt(u, p, t; model=model, st=st)
    model(u, p, st)[1]
end

prob = ODEProblem(dudt, u0, tspan, saveat=tsteps)
sense = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))

function predict(ps; prob=prob, odealg=odealg, sense=sense)
    solve(prob,
          odealg,
          p=ps,
          sensealg=sense
         ) |> CuArray
end

function loss(p)
    pred = predict(p)
    loss = sum(abs.(pred .- 1.0))

    loss, pred
end

function cb(p, l, pred; space=space)
    println(l)

    return false
end

# dummy
println("fwd"); cb(ps,loss(ps)...)
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
