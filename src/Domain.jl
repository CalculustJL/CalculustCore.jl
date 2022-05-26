#
###
# AbstractDomain interface #TODO
###

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

###
# IntervalDomain
###

""" 1D interval """
struct IntervalDomain{T<:Number,Ttag} <: AbstractDomain{T,1}
    x0::T
    x1::T
    periodic::Bool
    bdry_tags::Ttag

    function IntervalDomain(x0::Number, x1::Number, periodic, bdry_tags)
        T = promote_type(eltype.((x0, x1))...)
        new{T,typeof(bdry_tags)}(T(x0), T(x1), periodic, bdry_tags)
    end
end

function IntervalDomain(x0 = -1e0, x1 = 1e0;
                        periodic=false, bdry_tags = (nothing, nothing)
                       )
    IntervalDomain(x0, x1, periodic, bdry_tags)
end

function (::Type{T})(int::IntervalDomain) where{T<:Number}
    IntervalDomain(T(int.x0), T(int.x1), int.periodic, int.bdry_tags)
end

isperiodic(dom::IntervalDomain) = dom.periodic
endpoints(dom::IntervalDomain) = (dom.x0, dom.x1)
boundary_tags(dom::IntervalDomain) = dom.bdry_tags
boundary_tag(dom::IntervalDomain, i) = dom.bdry_tags[i]
num_boundaries(dom::IntervalDomain) = 2

function domains_match(int1::IntervalDomain, int2::IntervalDomain)
    bools = isapprox.(endpoints.((int1, int2))...)
    *(bools...)
end

###
# Box Domain
###

""" D-dimensional logically reectangular domain """
struct BoxDomain{T,D,Ti} <: AbstractDomain{T,D}
    intervals::Ti

    function BoxDomain(intervals::IntervalDomain...)
        T = promote_type(eltype.(intervals)...)
        D = length(intervals)
        new{T,D,typeof(intervals)}(intervals)
    end
end

function (::Type{T})(box::BoxDomain) where{T<:Number}
    BoxDomain(T.(box.intervals)...)
end

isperiodic(box::BoxDomain, dir::Integer) = box.intervals[dir] |> isperiodic
isperiodic(box::BoxDomain) = isperiodic.(box.intervals)

endpoints(box::BoxDomain, dir::Integer) = box.intervals[dir] |> endpoints
endpoints(box::BoxDomain) = endpoints.(box.intervals)

boundary_tags(box::BoxDomain, dir::Integer) = box.intervals[dir] |> boundary_tags
boundary_tags(box::BoxDomain) = boundary_tags.(box.intervals)

boundary_tag(box::BoxDomain, dir::Integer, i::Integer) = boundary_tags(box, dir)[i]
boundary_tag(box::BoxDomain, i) = boundary_tag(box, cld(i,2), 1 + rem(i-1, 2))

num_boundaries(box::BoxDomain{<:Number,D}) where{D} = 2D

Base.:*(int1::IntervalDomain, int2::IntervalDomain) = BoxDomain(int1, int2)
Base.:*(box1::BoxDomain, int2::IntervalDomain) = BoxDomain(box1.intervals..., int2)
Base.:*(int1::IntervalDomain, box2::BoxDomain) = BoxDomain(int1, box2.intervals...)
Base.:*(box1::BoxDomain, box2::BoxDomain) = BoxDomain(box1.intervals..., box2.intervals...)

function domains_match(box1::BoxDomain, box2::BoxDomain)
    if dims(box1) == dims(box2)
        @warn "dimension mismatch between $box1, $box2"
        return false
    end
    
    ret = true
    for i=1:D1
        ret * domains_match(box1.intervals[i], box2.intervals[i])
    end

    ret
end

###
# DeformedDomain
###

"""
Deform D-dimensional domain via mapping

args:
   -reference domain
   -mapping function
        (x1,...,xD) = map(r1, ..., rD)
   -isrescaling - do optimizations if
    mapping is a simple rescaling
        x1 = a1 + λ1 * x1(r1), ...,
        xD = aD + λD * xD(rD)
   -isseparable - do optimizations if
    mapping is separable
        x1 = x1(r1), ...,
        xD = xD(rD)
"""
struct DeformedDomain{T,D,Tdom<:AbstractDomain{T,D}, Tm} <: AbstractDomain{T,D}
    domain::Tdom
    mapping::Tm
    isseparable::Bool
    isrescaling::Bool
end

function deform(domain, mapping = nothing;
                isseparable = false, isrescaling = false
               )
    DeformedDomain(domain, mapping, isseparable, isrescaling)
end

function (::Type{T})(dom::DeformedDomain) where{T<:Number}
    BoxDomain(T(dom), dom.mapping, dom.isseparable, dom.isrescaling)
end
#
