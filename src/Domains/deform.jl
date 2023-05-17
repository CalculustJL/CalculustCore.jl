#
###
# DomainMap
###
"""
# Arguments
- reference domain

- mapping function
    (x1,...,xD) = map(r1, ..., rD)

- isrescaling - do optimizations if mapping is a simple rescaling
    x1 = a1 + λ1 × x1(r1),
    ...,
    xD = aD + λD × xD(rD)

-isseparable - do optimizations if mapping is separable
    x1 = x1(r1),
    ...,
    xD = xD(rD)

TODO - make types for AffineMap, LinearMap, SeparableMap, Transation, Rotation
TODO - write constructor that checks for method
"""
struct DomainMap{Tm}
    map::Tm
    isseparable::Bool
    isrescaling::Bool
end

function DomainMap(map; isseparable = false, isrescaling = false)
    DomainMap(map, isseparable, isrescaling)
end

###
# MappedDomain
###

"""
Deform D-dimensional domain via mapping
"""
struct MappedDomain{T, D, Tdom <: AbstractDomain{T,D}, Tmap<:DomainMap} <: AbstractDomain{T, D}
    domain::Tdom
    map::Tmap
end

function deform(domain::AbstractDomain, map = nothing;
                isseparable = false, isrescaling = false)
    isnothing(map) && return domain

    _map = DomainMap(map; isseparable = isseparable, isrescaling = isrescaling)

    MappedDomain(domain, _map)
end

function (::Type{T})(dom::MappedDomain) where {T <: Number}
    MappedDomain(T(dom.domain), mapping)
end

###
# Transfinite Interpolation
###

# function gordon_hall end #TODO

# """
#  Gordon Hall map - Transfinite interpolation
# """
# function mapBoxes(box1::BoxDomain{<:Number, D}, box2::BoxDomain{<:Number, D}) where {D}
#     mapping = nothing
# end
#
