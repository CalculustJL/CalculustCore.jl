#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, CUDA, Sundials

#N  = 4 
#f  = (u, p, t) -> zero(u)
#u0 = rand(N)
#ts = (0f0, 1f0)
#
##u0 = rand(N); alg = RadauIIA3() # OK
##u0 = rand(N) .+ im; alg = RadauIIA3() # Fails
#u0 = rand(N) .+ im; alg = RadauIIA3(autodiff=false) # OK
#u0 = rand(N) .+ im; alg = RadauIIA3(autodiff=false, linsolve=KrylovJL_GMRES())
#
##u0 = CUDA.rand(N); alg = RadauIIA3() # OK
##u0 = CUDA.rand(N); alg = CVODE_BDF() # Fails
#
#prob = ODEProblem(f, u0, ts)
#sol  = solve(prob, alg)
################################
###############################
N  = 4 
f1 = (u, p, t) -> -0.1 * u
f2 = (u, p, t) ->  0.0 * u
u0 = rand(N)
ts = (0f0, 1f0)

func = SplitFunction(f1, f2)

#u0 = rand(N); alg = RadauIIA3() # OK
u0 = rand(N) .+ im; alg = RadauIIA3(autodiff=false) # OK

prob = ODEProblem(func, u0, ts)
sol  = solve(prob, alg)
#
