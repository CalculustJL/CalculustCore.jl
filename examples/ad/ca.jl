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
using Lux, Random, SciMLSensitivity, Zygote

CUDA.allowscalar(false)

function dudt(u, p, t)
    zero(u)
end

N  = 128
u0 = ComponentArray(; vx=rand(N,10), vy=rand(N,10)) |> gpu
ps = ones(4) # unused

tspan  = (0f0, 1f0)
tsave  = range(tspan...; length=100)
prob   = ODEProblem(dudt, u0, tspan, ps; reltol=1f-4, abstol=1f-4)
sense  = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))
odealg = Tsit5()

function predict(p)
    sol = solve(prob, odealg, p=p, sensealg=sense, saveat=tsave)

    vxs = Tuple(sol.u[i].vx for i=1:length(sol))
    vx  = cat(vxs...;dims=3)

    vys = Tuple(sol.u[i].vy for i=1:length(sol))
    vy  = cat(vys...;dims=3)

    vx, vy
end

function loss(p)
    vx, vy = predict(p)
    loss = sum(abs2.(vx.- 1f0))

    loss, vx
end

# dummy calls
optf = p -> loss(p)[1]
println("fwd"); @time optf(ps) |> display                  # works
println("bwd"); @time Zygote.gradient(optf, ps) |> display # fails
#
