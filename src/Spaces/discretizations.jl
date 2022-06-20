#

"""
Weighted residual method
"""
struct Galerkin <: AbstractDiscretization end

_transp(a, ::Galerkin) = adjoint(a)

"""
Collocation
"""
struct Collocation <: AbstractDiscretization end

function massOp(space::AbstractSpace, discr::Collocation)
    N = length(space)
    IdentityOperator{N}()
end

_transp(a, ::Collocation) = _reshape(a, (1,length(a)))
#
