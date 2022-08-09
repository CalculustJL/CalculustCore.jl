#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, CUDA, SciMLOperators

N  = 4 
u0 = rand(N)
ts = (0f0, 1f0)

f1 = (u, p, t) -> -0.1 * u
f2 = (u, p, t) ->  0.0 * u
func = SplitFunction(f1, f2)

u0  = rand(N)
alg = RadauIIA3()

prob = ODEProblem(func, u0, ts)
sol  = solve(prob, alg)
#
#u0 = rand(N) .+ im; alg = RadauIIA3(autodiff=false) # OK
