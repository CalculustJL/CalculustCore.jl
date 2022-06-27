#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, LinearAlgebra
using Plots

N = 128
ν = 0e0
p = ()

""" space discr """
space = FourierSpace(N)
tspace = transform(space)
discr = Collocation()

(x,) = points(space)
(k,) = modes(space)
ftr  = transformOp(space)
itr  = transformOp(tspace)

""" operators """
Â = diffusionOp(ν, tspace, discr)

v = 1.0; vel = @. x*0 + v
Ĉ = advectionOp((ftr * vel,), tspace, discr)
F̂ = -Ĉ

Â = cache_operator(Â, k)
F̂ = cache_operator(F̂, k)

""" IC """
function uIC(x)
    @. sin(1x)
#   u0 = rand(ComplexF64, size(k))
#   u0[20:end] .= 0
#   ftr \ u0
end
u0 = uIC(x)
û0 = ftr * uIC(x)

""" time discr """
tspan = (0.0, 10.0)
tsave = (0, π/4, π/2, 3π/4, 2π,)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(Â, F̂, û0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
pred = itr.(sol.u, nothing, 0)
pred = hcat(pred...)

utrue(x,v,t) = uIC(@. x - v*t)
utr = utrue(x,v,sol.t[1])
for i=2:length(sol.u)
    ut = utrue(x,v,sol.t[i])
    global utr = hcat(utr, ut)
end

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=true)
    plot!(plt, x, utr[:,i], legend=true)
end
display(plt)

err = norm(pred .- utr,Inf)
display(err)
@test err < 1e-4
#
