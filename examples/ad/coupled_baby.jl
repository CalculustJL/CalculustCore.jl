#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearAlgebra
using SciMLSensitivity, Zygote
using Plots

N  = 2
Nt = 4

u0 = [1f0, 1f0]
ps = [-1f0]

A  = [-1f0 1f0
       0f0 0f0]

dudt = function(u, p, t)

    A * u + [0, p[1]]
end

tspan  = (0f0, 1f0)
tsave  = range(tspan...; length=Nt)
prob   = ODEProblem(dudt, u0, tspan, ps; reltol=1f-4, abstol=1f-4)
sense  = InterpolatingAdjoint(autojacvec=ZygoteVJP(allow_nothing=true))
odealg = Tsit5()

predict = function(p)
    solve(prob, odealg, p=p, sensealg=sense, saveat=tsave)
end

loss = function(p)
    pred = predict(p) |> Array
    vx = @views pred[1, :]
    vy = @views pred[2, :]

    loss = sum(abs2.(vx.- 0f5))

    loss, vx
end

# dummy calls
optf = p -> loss(p)[1]
println("fwd"); @time optf(ps) |> display                  # works
println("bwd"); @time Zygote.gradient(optf, ps) |> display # works
#
