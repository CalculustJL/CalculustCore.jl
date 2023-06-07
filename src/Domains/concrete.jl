#
function periodic_interval_tags(dim::Integer)
    tag = Symbol("D$(dim)_periodic")

    (tag, tag)
end

function default_interval_tags(dim::Integer)
    tag1 = Symbol("D$(dim)_inf")
    tag2 = Symbol("D$(dim)_sup")

    (tag1, tag2)
end

UnitBoxDomain(D; kwargs...) = FixedEndpointBoxDomain(D, 0.0, 1.0; kwargs...)
ChebyshevDomain(D; kwargs...) = FixedEndpointBoxDomain(D, -1.0, 1.0; kwargs...)
FourierDomain(D; kwargs...) = FixedEndpointBoxDomain(D, -π, π; periodic_dims = 1:D, kwargs...)

UnitIntervalDomain(; kwargs...) = UnitBoxDomain(1; kwargs...)
UnitSquareDomain(; kwargs...) = UnitBoxDomain(2; kwargs...)
UnitCubeDomain(; kwargs...) = UnitBoxDomain(3; kwargs...)

"""
`D`-dimensional `BoxDomain` with same endpoints for all dimensions.
"""
function FixedEndpointBoxDomain(D::Integer, x0::Number, x1::Number; kwargs...)
    endpoints = ()

    for _ in 1:D
        endpoints = (endpoints..., x0, x1)
    end

    BoxDomain(endpoints...; kwargs...)
end

"""
Generate a `D`-dimensional box domain with endpoints given by `endpoints`.
"""
function BoxDomain(endpoints::Real...;
                   periodic_dims = (),
                   tag = nothing,
                   boundary_tags = nothing)

    @assert length(endpoints) |> iseven "Number of endpoints must be even!"
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

IdentityMap() = DomainMap(identity_map; isseparable = true, isrescaling = true)
PolarMap() = DomainMap(polar_map, isseparable = false, isrescaling = false)

function AnnulusDomain(rIn, rOut)
    intR = IntervalDomain(rIn, rOut, false, (:Inner, :Outer))
    intθ = IntervalDomain(-π, π, true, (:Periodic, :Periodic))

    dom = intR × intθ

    deform(dom, PolarMap())
end
#
