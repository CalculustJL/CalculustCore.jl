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
check if domain end points match

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
    isperiodic::Bool
    tags::Ttag

    function IntervalDomain(x0::Number, x1::Number, isperiodic, tags)
        T = promote_type(eltype.((x0, x1))...)
        new{T,typeof(tags)}(T(x0), T(x1), isperiodic, tags)
    end
end

function IntervalDomain(x0 = -1e0, x1 = 1e0,
                        isperiodic=false, tags = (nothing, nothing)
                       )
    IntervalDomain(x0, x1, isperiodic, tags)
end

isperiodic(dom::IntervalDomain) = dom.isperiodic
endpoints(dom::IntervalDomain) = (dom.x0, dom.x1)

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

isperiodic(box::BoxDomain, dir::Integer) = box.intervals[dir] |> isperiodic
endpoints(box::BoxDomain, dir::Integer) = box.intervals[dir] |> endpoints
endpoints(box::BoxDomain) = endpoints.(box.intervals)

Base.:*(int1::IntervalDomain, int2::IntervalDomain) = BoxDomain(int1, int2)
Base.:*(box1::BoxDomain, int2::IntervalDomain) = BoxDomain(box1.intervals..., int2)
Base.:*(int1::IntervalDomain, box2::BoxDomain) = BoxDomain(int1, box2.intervals...)
Base.:*(box1::BoxDomain, box2::BoxDomain) = BoxDomain(box1.intervals..., box2.intervals...)

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
#   isrescaling::Bool
end

function deform(domain, mapping = nothing, isseparable = false)
    DeformedDomain(domain, mapping, isseparable)
end

function domains_match(box1::BoxDomain{<:Number, D1}, box2::BoxDomain{<:Number, D2}) where{D1,D2}
    if D1 != D2
        @warn "dimension mismatch between $box1, $box2"
        return false
    end
    
    ret = true
    for i=1:D1
        ret * domains_match(box1.intervals[i], box2.intervals[i])
    end

    ret
end

function map_from_ref(domain, ref_domain;D=D) # TODO

    if domains_match(domain, ref_domain)
        return domain
    end

    xe = domain.endpoints

    mapping = domain.mapping ==! nothing ? domain.mapping : (r -> r)
    mapping = mapping ∘ map
end
#
