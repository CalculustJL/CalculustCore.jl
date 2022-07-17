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
using Optimization, OptimizationOptimJL, OptimizationOptimisers
using Plots

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

""" fully oop problem """
#implicit(u, p, t) = D(u,p,t)
#function explicit(u, p, t; space=space, model=model, st=st)
#    x = points(space)[1]
#
#    dut = model(x', p, st)[1]
#    return vec(dut)
#end
#prob = SplitODEProblem{false}(implicit, explicit, u0, tspan, saveat=tsteps)
#sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())

""" fully iip problem """
#implicit(du, u, p, t) = D(du,u,p,t)
#function explicit(du, u, p, t; space=space, model=model, st=st)
#    x = points(space)[1]
#
#    dut = model(x', p, st)[1]
#    copyto!(du, dut)
#    return du
#end
#prob = SplitODEProblem{true}(implicit, explicit, u0, tspan, saveat=tsteps)
#sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())

""" total """
c1 = copy(x)
function ddt(du, u, p, t; space=space, model=model, st=st, c1=c1)
    x = points(space)[1]

    # implicit
    D(c1, u, p, t)

    # explicit
    dut = model(x', p, st)[1]
    copyto!(du, dut)

    axpy!(true, c1, du)

    return du
end
prob = ODEProblem{true}(ddt, u0, tspan, saveat=tsteps)
sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())

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
println("fwd"); cb(ps,loss(ps)...;doplot=false)
println("bwd"); Zygote.gradient(p -> loss(p)[1], ps) |> display

#""" optimization """
#adtype = Optimization.AutoZygote()
#optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype) # x=object to optimize
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
