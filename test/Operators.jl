#
using PDEInterfaces, LinearAlgebra

# OperatorBasics.jl
import PDEInterfaces: NullOp, IdentityOp, AffineOp, ComposedOp, InverseOp

# Operators.jl
import PDEInterfaces: MatrixOp, DiagonalOp, TensorProductOp2D

nr = 8
ns = 12

p = nothing
t = false

u = rand(nr,ns) |> Field
v = rand(nr,ns) |> Field
w = rand(nr,ns) |> Field

# DiagonalOp
d = rand(nr,ns) |> Field
D = DiagonalOp(d)
D_u = d .* u

# TensorProductOp2D
A = rand(nr,nr)
B = rand(ns,ns)
T = TensorProductOp2D(A, B)
T_u = Field(A * u.array * B')

# NullOp
Z = NullOp{2}()
Z_u = u * false

# IdentityOp
Id = IdentityOp{2}()
Id_u = copy(u)

ops = (
       Z, Id, D, T,
      )

rhs = (
       Z_u, Id_u, D_u, T_u,
      )

for (A, b) in zip(ops, rhs)
    # in place
    @test mul!(v, A, u) ≈ b
    @test A(v, u, p, t) ≈ b

    # out of place
    @test *(A, u) ≈ b
#   @test A(u, p, t) ≈ b
end

## AffineOp

# +,-(op,op)
# +,-(λ, op)
# *,/(λ, op)

## ComposeOp

# #TODO caching mechanism is broken. fix later. 

## InverseOp

