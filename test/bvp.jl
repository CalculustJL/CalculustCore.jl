#
using PDEInterfaces
using LinearSolve

domain = reference_box(1)
space = GaussLobattoLegendre1D(32; domain=domain)
(x,) = grid = get_grid(space)

op = laplaceOp(space)
f  = @. 0*x + 1
bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => NeumannBC(),
          )

prob = BoundaryValuePDEProblem(op, f, bcs, space)
alg  = LinearBVPDEAlg()
sol  = solve(prob, alg)
