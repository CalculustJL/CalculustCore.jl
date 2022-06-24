#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, Sundials, LinearAlgebra
using Zygote, Random, Lux, DiffEqSensitivity, ComponentArrays

N = 128
ν = 1e-1

odealg = Tsit5()
#odealg = CVODE_BDF(method=:Functional)
#odealg = Rodas5(autodiff=false)

""" space discr """
domain = FourierDomain(1)
space  = FourierSpace(N; domain=domain)
discr  = Collocation()

(x,) = points(space)

# diffusion
D = diffusionOp(ν, space, discr)

# forcing
Z = SciMLOperators.NullOperator{length(space)}()
F = Z
#F = AffineOperator(Z, zero(x)) <-- NN output

D = cache_operator(D, x)
F = cache_operator(F, x)

""" NN """
model = Lux.Chain(
                  Lux.Dense(1, 10),
                  Lux.Dense(10, 1),
                 )

rng = Random.default_rng()
ps, st = Lux.setup(rng, model)
ps = ComponentArray(ps)

u0 = @. sin(10x)
tspan = (0.0, 1.0)
tsteps = range(tspan..., length=10)

""" fully explicit problem """
function implicit(u, p, t;op=D)
    op * u
end

function explicit(u, p, t; space=space, model=model, st=st)
    x  = points(space)[1]
    xt = transpose(x)

    dut = model(xt, p, st)[1]
    return vec(dut)
end

#prob = SplitODEProblem{false}(implicit, explicit, u0, tspan, saveat=tsteps)
#sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())

prob = ODEProblem{false}(explicit, u0, tspan, saveat=tsteps)
sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())

function predict(ps; prob=prob, odealg=odealg, sense=sense)
    solve(prob,
          odealg,
          p=ps,
          sensealg=sense
         ) |> Array
end

function loss(p)
    pred = predict(p)
    loss = sum(abs.(pred .- 1.0))

    loss, pred
end

function cb(p, l, pred; doplot=true, space=space)
    println(l)
    return false
end

# dummy
println("fwd"); cb(ps,loss(ps)...;doplot=true)
println("bwd"); Zygote.gradient(p -> loss(p)[1], ps)

#=
""" optimization """
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss_neuralode(x), adtype)
optprob = Optimization.OptimizationProblem(optf, pinit)

optres = Optimization.solve(optprob,
                            ADAM(0.05),
                            callback = callback,
                            maxiters = 300
                           )

optprob = remake(optprob,u0 = result_neuralode.u)

optres = Optimization.solve(optprob,
                            Optim.BFGS(initial_stepnorm=0.01),
                            callback=callback,
                            allow_f_increases = false
                           )

callback(optres.u, loss(optres.u)...; doplot=true)
=#
#
