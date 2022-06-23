#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces
using OrdinaryDiffEq, LinearSolve, Sundials, LinearAlgebra
using Zygote, Random, Lux, DiffEqSensitivity

N = 128
ν = 1e-2
p = ()
odealg = Tsit5()
odealg = CVODE_BDF(method=:Functional)
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
rng = Random.default_rng()

model = Lux.Chain(
                  Lux.Dense(1, 10),
                  Lux.Dense(10, 1),
                 )

ps, st = Lux.setup(rng, model)
ps = ComponentArray(ps)

function explicit!(du, u, p, t; space=space, model=model, st=st)
    x  = points(space)[1]
    xt = transpose(x)

    dut = model(xt, p, st)[1]
    du .= transpose(dut)
end

""" IC """
u0 = @. sin(10x)

""" solve """
tspan = (0.0, 10.0)
tsteps = range(tspan..., length=100)

prob = SplitODEProblem{true}(D, explicit!, u0, tspan, ps, saveat=tsteps)

sense = InterpolatingAdjoint(autojacvec=ZygoteVJP())

function loss(p; prob=prob, odealg=odealg)
    prob = remake(prob, p=p)
    pred = solve(prob, odealg) |> Array
    loss = sum(abs.(pred .- 1.0))
    loss, pred
end

function cb(p, l, pred; doplot=true, space=space)
    println(l)
#   if doplot
#       x = points(space)[1]
#       plot(x, )
#   end
    return false
end

# dummy
println("fwd"); cb(ps,loss(ps)...;doplot=true)
println("bwd"); Zygote.gradient(p -> loss(p)[1], p)

#=
""" fully explicit problem """
function implicit(u, p, t;op=D)
    op * u
end

function explicit(u, p, t; space=space, model=model, st=st)
    x  = points(space)[1]
    xt = transpose(x)

    dut = model(xt, p, st)[1]
    return transpose(dut)
end

prob = SplitODEProblem{false}(implicit, explicit, u0, tspan, ps, saveat=tsteps)
=#

#
nothing
