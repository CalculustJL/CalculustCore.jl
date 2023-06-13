#
###
# Collocation
###

COLLOCATION_DEALIAS_MSG = """
    When forming nonlinear operators...
    """

"""
$TYPEDEF

A `Collocation` discretization is when the strong form of the equation is
    discretized as it is.
"""
struct Collocation <: AbstractDiscretization end

"""
$SIGNATURES

```
[Dx  --> [Dx Dy]
 Dy]
```
"""
_transp(v::AbstractVector, ::Collocation) = reshape(v, (1, length(v)))

"""
For a `Collocation` discretization, the mass operator returns
`IdentityOperator(V)`.
"""
massOp(V::AbstractSpace, ::Collocation) = IdentityOperator(V)

function laplaceOp(V::AbstractSpace, ::Collocation)
    DD2 = hessianOp(V)
    -tr(DD2) # trace
end

function biharmonicOp(V::AbstractSpace, discr::Collocation)
    A = laplaceOp(V, discr)

    A * A
end

"""
$SIGNATURES

$COLLOCATION_DEALIAS_MSG
"""
function advectionOp(vels::NTuple{D},
                     V::AbstractSpace{<:Any, D},
                     discr::Collocation,
                     Vd::AbstractSpace{<:Any, D},
                     J = nothing;
                     vel_update_funcs = nothing) where {D}

    @error "TODO: dealiased collocation advection has not been implemented"

    # J12 = J !== nothing ? J : interpOp(V1, V2)
    #
    # VV1 = [DiagonalOperator.(vels)...]
    # VV2 = J12 .* VV1
    # M2 = massOp(V2, discr)
    # MM2 = Diagonal([M for i in 1:D])
    # DD1 = gradientOp(V1)
    #
    # VV2' * MM2 * (J12 .* DD1)
end

function advectionOp(vels::NTuple{D},
                     V::AbstractSpace{<:Any, D},
                     discr::Collocation,
                     truncation_fracs::NTuple{D};
                     vel_update_funcs = nothing) where {D}

    @error "TODO: collocation advection with spectral truncation has not
    been implemented"

    # TODO: truncation same as dealiasing to a lower order space

    # generalize `FourierSpaces.jl` implementation
end

function divergenceOp(V::AbstractSpace, ::Collocation)
    DD = gradientOp(V)

    sum(DD)
end

