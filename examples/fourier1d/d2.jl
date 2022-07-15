using PDEInterfaces, OrdinaryDiffEq, Plots
N = 32

domain = GLLBox(2)
space = GaussLobattoLegendre(N, N; domain=domain)
(x, y,) = points(space)

op = laplaceOp(space)
f  = @. 0*x + 1
bcs = Dict(
           :Lower1 => NeumannBC(),
           :Upper1 => DirichletBC(),

           :Lower2 => DirichletBC(),
           :Upper2 => NeumannBC(),
          )

prob = BVPDEProblem(op, f, bcs, space)
alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
plt = plot(sol)
savefig(plt, "bvp2d_dn")
