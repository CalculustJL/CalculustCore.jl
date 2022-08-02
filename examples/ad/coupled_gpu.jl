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

N  = 2
K  = 3
Nt = 4

u0 = ComponentArray(; vx=ones(Float32, N, K), vy=ones(Float32, N, K)) |> gpu
ps = [1f0] |> gpu

dudt = function(u, p, t)
    dvx = -1f0 * u.vx +  1f0 * u.vy
    dvy = p   .* u.vy

    dvx = Zygote.hook(Δ -> (println("||Δdvx||∞ ", norm(Δ, Inf)); Δ), dvx)
    dvy = Zygote.hook(Δ -> (println("||Δdvy||∞ ", norm(Δ, Inf)); Δ), dvy)

    ComponentArray(vcat(dvx |> vec, dvy |> vec), getaxes(u))
end

tspan  = (0f0, 1f0)
tsave  = range(tspan...; length=Nt)
prob   = ODEProblem(dudt, u0, tspan, ps; reltol=1f-4, abstol=1f-4)
sense  = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))
odealg = Tsit5()

predict = function(p)
    pred = solve(prob, odealg, p=p, sensealg=sense, saveat=tsave) |> CuArray

    pred = Zygote.hook(Δ -> (println("Δpred isa ", typeof(Δ)); Δ), pred)

    vx = @views pred[1    : N*K, :] # size [N,K,Nt]
    vy = @views pred[N*K+1:2N*K, :]

    vx = reshape(vx, (N,K,Nt))
    vy = reshape(vy, (N,K,Nt))

    #vx = Zygote.hook(Δ -> (println("Δvx isa ", typeof(Δ)); Δ), vx)
    #vy = Zygote.hook(Δ -> (println("Δvy isa ", typeof(Δ)); Δ), vy)

    vx, vy
end

loss = function(p)
    vx, vy = predict(p)
    loss = sum(abs2.(vx.- 0f5))

    loss, vx
end

# dummy calls
optf = p -> loss(p)[1]
println("fwd"); @time optf(ps) |> display                  # works
println("bwd"); @time Zygote.gradient(optf, ps) |> display # works
gr = Zygote.gradient(optf, ps)[1]
#
