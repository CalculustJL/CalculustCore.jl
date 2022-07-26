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
using DiffEqSensitivity, Zygote

N = 128
ν = 1e-1
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
domain = FourierDomain(1)
space  = FourierSpace(N; domain=domain)
space  = make_transform(space; p=ps)
discr  = Collocation()

(x,) = points(space)

D = diffusionOp(ν, space, discr)
D = cache_operator(D, x)

u0 = @. sin(10x)
tspan = (0.0, 1.0)
tsteps = range(tspan..., length=10)

""" fully iip problem """
implicit(du, u, p, t) = D(du,u,p,t)
function explicit(du, u, p, t; space=space, model=model, st=st)
    x = points(space)[1]

    dut = model(x', p, st)[1]
    copyto!(du, dut)
    return du
end

function vjp(Jv,v,u,p,t)
    return Jv
end

odefunc = SplitFunction{true}(implicit, explicit; vjp=vjp)
odefunc = SplitFunction{true}(implicit, explicit; vjp=vjp)

prob = ODEProblem(odefunc, u0, tspan, saveat=tsteps)
sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())
#sense = InterpolatingAdjoint(autojacvec=ReverseDiffVJP())

function predict(ps; prob=prob, odealg=odealg, sense=sense)
    solve(prob, odealg, p=ps, sensealg=sense) |> Array
end

function loss(p)
    pred = predict(p)
    loss = sum(abs.(pred .- 1.0))

    loss, pred
end

# dummy
println("fwd"); loss(ps)[1] |> display
println("bwd"); Zygote.gradient(p -> loss(p)[1], ps) |> display
#
