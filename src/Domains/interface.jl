#
Base.eltype(::AbstractDomain{T}) where {T} = T

#Base.in
#Base.:∈

"""
    dims(dom)

Dimension of domain
"""
dims(::AbstractDomain{<:Any, D}) where {D} = D

"""
    bounding_box(dom)

Returns bounding box of domain as a BoxDomain of the same dimension.

Must overload.
"""
function bounding_box end

"""
    expanse(dom)

Returns size of domain bounding box.

Overload if `expanse` is not defined for `bounding_box(dom)`.
"""
expanse(dom::AbstractDomain) = expanse(bounding_box(dom))

"""
    isperiodic(dom)

Returns NTuple{D, Bool} indicating if the domain in periodic along each
dimension.
"""
function isperiodic end

"""
    boundaries(dom)

Returns tuple of domain boundaries of `dom` as `D-1` domains.

Must overload.
"""
function boundaries end

"""
    domain_tag(dom)

Returns metadata tag corresponding to the interior of domain `dom`.

Must overload
"""
function domain_tag end

"""
num_boundaries(dom)

Returns number of boundaries in domain `dom`.
"""
num_boundaries(dom::AbstractDomain) = length(boundaries(dom))

"""
    boundary_tags(dom)

Returns tuple of tags corresponding to `boundaries(dom)`.
"""
function boundary_tags(dom::AbstractDomain)
    domain_tag.(boundaries(dom))
end
#
