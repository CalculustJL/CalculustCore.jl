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
    expanse(dom)

Returns NTuple{D} of domain lengths along each dimension
"""
expanse(dom::AbstractDomain) = expanse(bounding_box(dom))

"""
    isperiodic(dom)

Returns NTuple{D, Bool} indicating if the domain in periodic along each
dimension.
"""
function isperiodic(dom::AbstractDomain{T, D}) where{T, D}
    Tuple(isperiodic(dom, d) for d in 1:D)
end

"""
    isperiodic(dom, dim)

Returns Bool indicating if `dom` is periodic in dimension `d`.
"""
function isperiodic end

"""
    boundaries(dom)

Returns tuple of domain boundaries of `dom` as a `D-1` domains.
"""
function boundaries end

# """
#     boundary(dom)
#
# Returns boundary of domain `dom` as a `D-1` domain.
# """
# boundary(dom::AbstractDomain) = ∪(boundaries(dom)...)

"""
num_boundaries(dom)

Returns number of boundaries in domain `dom`.
"""
num_boundaries(dom::AbstractDomain) = length(boundaries(dom))

"""
    domain_tag(dom)

Returns metadata tag corresponding to the interior of domain `dom`.
"""
function domain_tag end

"""
    boundary_tag(dom)

Returns tuple of tags corresponding to `boundaries(dom)`.
"""
function boundary_tag(dom::AbstractDomain)
    domain_tag.(boundaries(dom))
end

"""
    boundary_tag(dom, i)

returns tuple of metadata tags corresponding to `boundaries(dom)`.
"""
function boundary_tag end

"""
    bounding_box(dom)

Returns bounding box of domain as a BoxDomain of the same dimension.
"""
function bounding_box end
#
