#
"""
$TYPEDEF

Logically rectangular product domain
"""
const BoxedDomain{T, D} = ProductDomain{T, D, <:NTuple{D, IntervalDomain}} where{T, D}

const LogicallyRectangularDomain{T,D} = Union{
                                              IntervalDomain{T},
                                              BoxedDomain{T, D},
                                              MappedDomain{T, D, <:BoxedDomain},
                                             } where{T, D}

function periodic_interval_tags(dim::Integer)
    tag = Symbol("D$(dim)_periodic")

    (tag, tag)
end

function default_interval_tags(dim::Integer)
    tag1 = Symbol("D$(dim)_inf")
    tag2 = Symbol("D$(dim)_sup")

    (tag1, tag2)
end

BOX_KW_MSG = """
    # Keyword Arguments
    - `periodic_dims`: List of periodic dimensions in `1:D`. Defaults to none.
    - `tag`: Domain tag (symbol) defaults to `:Interior`.
    - `boundary_tags`: List of the `2D` boundary tags. By default, periodic boundaries are tagged `:Dk_periodic` where `k` refers to the `k`th dimension, and aperiodic boundaries are tagged as `:Dk_inf`, `:Dk_sup` corresponding to the infimum and supremum in the `k`th dimension.
    """

"""
$SIGNATURES

Creates a `D`-dimensional, unit sized box domain with endpoints at `0.0`, `1.0`
in each dimension.

$BOX_KW_MSG
"""
UnitBoxDomain(D; kwargs...) = FixedEndpointBoxDomain(D, 0.0, 1.0; kwargs...)

"""
$SIGNATURES

Creates a `D`-dimensional box domain with endpoints at `-1.0`, `1.0` in each
dimension.

$BOX_KW_MSG
"""
ChebyshevDomain(D; kwargs...) = FixedEndpointBoxDomain(D, -1.0, 1.0; kwargs...)

"""
$SIGNATURES

Creates a `D`-dimensional box domain that is periodic in all directions.
The domain spans from `-π` to `π` in each dimension.

$BOX_KW_MSG
"""
FourierDomain(D; kwargs...) = FixedEndpointBoxDomain(D, -π, π; periodic_dims = 1:D, kwargs...)

"""
$SIGNATURES

Creates a unit-sized 1D interval domain spanning from `0.0` to `1.0`.

$BOX_KW_MSG
"""
UnitIntervalDomain(; kwargs...) = UnitBoxDomain(1; kwargs...)

"""
$SIGNATURES

Creates a unit-sized 2D square domain spanning from `0.0` to `1.0` in each
dimension.

$BOX_KW_MSG
"""
UnitSquareDomain(; kwargs...) = UnitBoxDomain(2; kwargs...)

"""
$SIGNATURES

Creates a unit-sized 3D cube domain spanning from `0.0` to `1.0` in each
dimension.

$BOX_KW_MSG
"""
UnitCubeDomain(; kwargs...) = UnitBoxDomain(3; kwargs...)

"""
$SIGNATURES

Creates a `D`-dimensional box domain with the same end points in each dimension.

$BOX_KW_MSG
"""
function FixedEndpointBoxDomain(D::Integer, x0::Number, x1::Number; kwargs...)
    endpoints = ()

    for _ in 1:D
        endpoints = (endpoints..., x0, x1)
    end

    BoxDomain(endpoints...; kwargs...)
end

"""
$SIGNATURES

Creates a `D`-dimensional box domain with endpoints given by `endpoints`

$BOX_KW_MSG

# Example

```
julia> domain = BoxDomain(0, 1, 2, 3; periodic_dims=2)
(0, 1) × (2, 3)

julia> boundary_tags(domain)
(:D1_inf, :D1_sup, :D2_periodic, :D2_periodic)
```
"""
function BoxDomain(endpoints::Real...;
                   periodic_dims = (),
                   tag = nothing,
                   boundary_tags = nothing)

    if !iseven(length(endpoints))
        msg = """Number of endpoints must be even!"""
        throw(ArgumentError(msg))
    end

    D = div(length(endpoints), 2)
    tag = isnothing(tag) ? :Interior : tag

    intervals = ()
    for d in 1:D
        idx = [2d - 1, 2d]

        periodic = d ∈ periodic_dims
        bdr_tags = if isnothing(boundary_tags)
            periodic ? periodic_interval_tags(d) : default_interval_tags(d)
        else
            boundary_tags[idx]
        end

        interval = IntervalDomain(endpoints[idx]...;
                                  periodic = periodic,
                                  boundary_tags = bdr_tags)

        intervals = (intervals..., interval)
    end

    ProductDomain(intervals...; tag = tag)
end

###
# Maps
###

function polar_map(r::AbstractArray, θ::AbstractArray)
    x = @. r * cos(θ)
    y = @. r * sin(θ)
    x, y
end

identity_map(xyz::AbstractArray...) = xyz

"""
$SIGNATURES

The Identity `DomainMap`
"""
IdentityMap() = DomainMap(identity_map; isseparable = true, isrescaling = true)

"""
$SIGNATURES

2D Polar transformation
"""
PolarMap() = DomainMap(polar_map, isseparable = false, isrescaling = false)

"""
$SIGNATURES

Annulus domain with internal radius `r0`, and external radius `r1`.
"""
function AnnulusDomain(r0, r1)
    intR = IntervalDomain(r0, r1, false, (:Inner, :Outer))
    intθ = IntervalDomain(-π, π, true, (:Periodic, :Periodic))

    dom = intR × intθ

    deform(dom, PolarMap())
end
#
