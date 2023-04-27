#
using PDEInterfaces
using LinearAlgebra, LinearSolve, SciMLOperators

N = 32
dom = FourierDomain(1)
discr = Collocation()

space = FourierSpace(N)

(x,) = pts = points(space)
D = gradientOp(space, discr) |> first

u0 = @. 0 * x + 1
u1 = @. sin(10 * x)

# derivative
D = cache_operator(D, x)
v = copy(x)
@test ≈(mul!(v, D, u0), @. 0.0 * x; atol = 1e-8)
@test ≈(mul!(v, D, u1), @. (10) * cos(10 * x); atol = 1e-8)

# interpolation

# integration

#
