#
include("NDgrid.jl")
using SciMLOperators: ⊗, IdentityOperator, _reshape, _vec

_transp(a, ::AbstractDiscretization) = transpose(a)
