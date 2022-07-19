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
using Plots, Test

N = 128
ν = 0e0
p = nothing

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

Â = cache_operator(Â, im*k)
F̂ = cache_operator(F̂, im*k)

""" IC """
X  = truncationOp(space, (1//2,))
u0 = X * rand(size(x)...)
û0 = ftr * u0

""" time discr """
tspan = (0.0, 10.0)
tsave = (0, π/4, π/2, 3π/4, 2π,)
odealg = Rodas5(autodiff=false)
prob = SplitODEProblem(Â, F̂, û0, tspan, p)
@time sol = solve(prob, odealg, saveat=tsave)

""" analysis """
pred = itr.(sol.u, nothing, 0)
pred = hcat(pred...)

#utrue(x,v,t) = uIC(@. x - v*t)
#utr = utrue(x,v,sol.t[1])
#for i=2:length(sol.u)
#    ut = utrue(x,v,sol.t[i])
#    global utr = hcat(utr, ut)
#end

plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, pred[:,i], legend=true)
end
display(plt)

#err = norm(pred .- utr,Inf)
#display(err)
#@test err < 1e-4
#
