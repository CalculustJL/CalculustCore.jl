#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using OrdinaryDiffEq, DiffEqSensitivity, LinearAlgebra
using Flux, Lux, Random, ComponentArrays

## FLUX
model = Flux.Chain(
                  Flux.Dense(4, 10),
                  Flux.Dense(10, 4),
                 )

p, re = Flux.destructure(model)
dudt(u, p, t; model=model, re=re) = re(p)(u)

## LUX
#model = Lux.Chain(
#                  Lux.Dense(4, 10),
#                  Lux.Dense(10, 4),
#                 )
#
#rng = Random.default_rng()
#p, st = Lux.setup(rng, model)
#dudt(u, p, t; model=model, st=st) = model(u, p, st)[1]

###

p = ComponentArray(p)

u0 = rand(4)
tsp = (0.0, 1.0)
tsv = Array(0:0.5:1.0)
prob = ODEProblem(dudt, u0, tsp)

function predict(p)
    solve(prob,
          Tsit5(),
          saveat=tsv,
          p=p,
#         sensealg=InterpolatingAdjoint(autojacvec=ZygoteVJP()),
          sensealg=InterpolatingAdjoint(autojacvec=EnzymeVJP()),
         ) |> Array
end

function loss(p)
    pred = predict(p)
    loss = sum(abs2, pred .- 1.0)

    loss, pred
end

println("fwd"); loss(p)[1] |> display
println("bwd"); gradient(p -> loss(p)[1], p) |> display
