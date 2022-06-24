#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using SciMLOperators, LinearSolve

N = 32

dom = reference_box(1)
space = GaussLobattoLegendre1D(N; domain=dom)
(x,) = pts = grid(space)

op = laplaceOp(space)
f  = @. 0*x + 1
bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => DirichletBC(),
          )

prob = BVPDEProblem(op, f, bcs, space)
alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
plt = plot(sol)
savefig(plt, "bvp_dd")
