#
using PDEInterfaces, LinearAlgebra

domain = reference_box(1)
space = GaussLobattoLegendre1D(32; domain=domain)

bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => NeumannBC(),
          )

#discr = GalerkinDiscretization()

(x,) = get_grid(space)
lhsOp = laplaceOp(space)
f = @. x * x

#
domain = reference_box(1)
space = GaussLobattoLegendre1D(32; domain=domain)

bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => NeumannBC(),
           :Lower2 => DirichletBC(),
           :Upper2 => NeumannBC(),
          )

(x,) = get_grid(space)
lhsOp = laplaceOp(space)
f = @. x * x

#prob = BoundaryValueProblem()

