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
using SciMLSensitivity, Zygote, Lux

CUDA.allowscalar(false)

function dudt(u, p, t)
    zero(u)
end

N  = 2
K  = 3
Nt = 4

u0 = ComponentArray(; vx=ones(N,K), vy=zeros(N,K)) |> gpu
ps = ones(1)

tspan  = (0f0, 1f0)
tsave  = range(tspan...; length=Nt)
prob   = ODEProblem(dudt, u0, tspan, ps; reltol=1f-4, abstol=1f-4)
sense  = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))
odealg = Tsit5()

predict = function(p)
    sol  = solve(prob, odealg, p=p, sensealg=sense, saveat=tsave)
    pred = sol |> CuArray

    NK = N*K
    vx = @views pred[1   : NK, :] # size [N,K,Nt]
    vy = @views pred[NK+1:2NK, :]

    vx, vy
end

loss = function(p)
    vx, vy = predict(p)
    loss = sum(abs2.(vx.- 1f0))

    loss, vx
end

# dummy calls
optf = p -> loss(p)[1]
println("fwd"); @time optf(ps) |> display                  # works
println("bwd"); @time Zygote.gradient(optf, ps) |> display # works
#
