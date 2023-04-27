#

Base.eltype(::AbstractDomain{T}) where {T} = T

#Base.in
#Base.:∈

"""
Dimension of domain
"""
dims(::AbstractDomain{<:Any, D}) where {D} = D

"""
Length of domain
"""
lengths(dom::AbstractDomain) = lengths(bounding_box(dom))

"""
args:
    AbstractDomain{T,D}
    direction < D
ret:
    Bool
"""
function isperiodic end

"""
args:
    AbstractDomain
    direction
ret:
    Tuple of end points
"""
function endpoints end

"""
args:
    AbstractDomain
    direction
ret:
    Tuple of boundary tags
"""
function boundary_tags end

"""
args:
    AbstractDomain
    i
ret:
    Tag of ith boundary
"""
function boundary_tag end

"""
get number of boundaries

args:
    AbstractDomain
ret:
    Integer
"""
function num_boundaries end

"""
get bounding box for domain

args:
    AbstractDomain
ret:
    BoxDomain
"""
function bounding_box end

"""
check if domain extent matches.
doesn't check periodicity or bdry_tags

args:
    AbstractDomain
    direction
ret:
    Bool
"""
function domains_match end
