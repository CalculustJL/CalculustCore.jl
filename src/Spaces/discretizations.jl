#
"""
$TYPEDEF

The `Galerkin` discretization is a weighted residual method where the trial
space is equal to the test space.
"""
struct Galerkin <: AbstractDiscretization end

_transp(a, ::Galerkin) = adjoint(a)

function laplaceOp(space::AbstractSpace, discr::Galerkin)
    D = ndims(space)

    M = massOp(space, discr)
    MM = Diagonal([M for i in 1:D])

    DD = gradientOp(space, discr)

    -DD' * MM * DD
end

"""
for v,u in H¹₀(Ω)

(v,-∇² u) = (vx,ux) + (vy,uy)\n
         := a(v,u)\n
"""
function laplaceOp(V1::AbstractSpace{<:Any, D},
                   V2::AbstractSpace{<:Any, D},
                   discr::Galerkin;
                   J = nothing) where {D}

    J12 = J !== nothing ? J : interpOp(V1, V2)
    #J21 = _transp(J12) # or interpOp(V2, V1) # TODO

    M2 = massOp(V2, discr)
    MM2 = Diagonal([M2 for i in 1:D])

    DD1 = gradientOp(V1, discr)
    JDD = J12 .* DD1

    JDD' * MM2 * JDD
end

function diffusionOp(ν::AbstractVector, space::AbstractSpace, discr::Galerkin)
    D = ndims(space)
    ν = DiagonalOperator(ν)
    DD = gradientOp(space)
    M = massOp(space, discr)
    Mν = ν * M
    MMν = Diagonal([Mν for i in 1:D])

    DD' * MMν * DD
end

"""
$TYPEDEF

A `Collocation` discretization
"""
struct Collocation <: AbstractDiscretization end
_transp(a, ::Collocation) = reshape(a, (1, length(a)))

"""
For a `Collocation` discretization, the mass operator returns
`IdentityOperator(V)`.

"""
massOp(V::AbstractSpace, ::Collocation) = IdentityOperator(V)

function laplaceOp(space::AbstractSpace, discr::Collocation)
    DD2 = -hessianOp(space, discr)
    sum(DD2)
end

function biharmonicOp(space::AbstractSpace, discr::Collocation)
    A = laplaceOp(space, discr)

    A * A
end

function divergenceOp(space::AbstractSpace, ::Collocation)
    DD = gradientOp(space)

    sum(DD)
end
#
