#
abstract type AbstractDiscretization end

struct Galerkin <: AbstractDiscretization end
struct Collocation <: AbstractDiscretization end

function massOp(space::AbstractSpace, discr::Collocation)
    N = length(space)
    IdentityOperator{N}()
end

function laplaceOp(space::AbstractSpace, discr::Collocation)
    D = gradOp(space, discr)

    D' * D
end

#
