#
using PDEInterfaces

domain = reference_box(1)
space = GaussLobattoLegendre1D(32; domain=domain)

op = laplaceOp(space)
f  = x -> @. 1 + x*0
bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => NeumannBC(),
          )

