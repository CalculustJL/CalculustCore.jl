#
using PDEInterfaces, LinearAlgebra

space = GaussLobattoLegendre1D(32)

bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => NeumannBC(),
          )

(x,) = get_grid(space)
lhsOp = laplaceOp(space)
f = @. x * x

#prob = BoundaryValueProblem()

