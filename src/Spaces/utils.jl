#
include("NDgrid.jl")
using SciMLOperators: ⊗, IdentityOperator, _reshape, _vec

_transp(a, ::AbstractDiscretization) = a'
_transp(a, ::Collocation) = _reshape(a, (1,length(a)))
