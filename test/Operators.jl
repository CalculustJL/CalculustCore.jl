#
using AbstractPDEs, LinearAlgebra

nr = 8
ns = 12

u = rand(nr,ns) |> Field
v = rand(nr,ns) |> Field

## DiagonalOp
d = rand(nr,ns) |> Field
D = DiagonalOp(d)

@test mul!(v,D,u) ≈ d .* u

## TensorProductOp2D

A = rand(nr,nr)
B = rand(ns,ns)
T = TensorProductOp2D(A, B)

@test mul!(v,T,u) ≈ Field(A * u.array * B')
@test mul!(v,T,u) ≈ Field(A * u.array * B')

## ZeroOp

## IdentityOp

## AffineOp

#### +,-(op,op)
#### +,-(λ, op)
#### *,/(λ, op)

## ComposeOp

# caching mechanism is broken. fix later. 

## InverseOp

