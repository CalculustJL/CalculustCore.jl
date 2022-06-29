#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, Lux
using Zygote, Random, DiffEqSensitivity, ComponentArrays

N = 10
w = 5
model = Lux.Chain(
                  Lux.Dense(N, w, tanh),
                  Lux.Dense(w, N),
                 )

rng = Random.default_rng()
ps, st = Lux.setup(rng, model)
ps = ComponentArray(ps)

function dudt(u, p, t; model=model, st=st)
    u_, st = model(u, p, st)
    u_
end

u0 = rand(N)
tspan = (0.0, 1.0)
tsave = range(tspan...; length=10)
prob = ODEProblem(dudt, u0, tspan, saveat=tsave)

function predict(p; prob=prob)
    solve(prob, Tsit5(), p=p) |> Array
end

function loss(p)
    pred = predict(p)
    loss = sum(abs2, pred .- 1.0)
    loss, pred
end

function callback(p, l, pred; doplot=false)
    println(l)
end

println("fwd"); callback(ps,loss(ps)...;doplot=true)
println("bwd"); Zygote.gradient(p -> loss(p)[1], ps)

#
