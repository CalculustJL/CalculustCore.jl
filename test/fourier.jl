#
using PDEInterfaces
using LinearAlgebra, LinearSolve, SciMLOperators

N = 32
dom = FourierDomain(1)
discr = Collocation()

space = FourierSpace(N)

(x,) = pts = points(space)
D = gradOp(space, discr) |> first

u0 = @. 0*x + 1
u1 = @. 1.0*x
u2 = @. sin(10 * x)
u4 = @. exp(x/π) / cos(x/π)
u5 = @. 1 / (1+16x^2)

# derivative
D = cache_operator(D, u)
v = copy(u)
@test ≈(mul!(v, D, u0), @. 0.0*x; atol=1e-8)
@test ≈(mul!(v, D, u1), @. 0.0*x + 1.0; atol=1e-8)
@test ≈(mul!(v, D, u2), @. (10) * cos(10 * x); atol=1e-8)

# interpolation

# integration

#
