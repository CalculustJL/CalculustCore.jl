#
###
# NullDomain
###

"""
Domain representing the empty set, ∅. Dimension set to `-1`.
"""
struct NullDomain <: AbstractDomain{Bool, -1} end
const ∅ = NullDomain()

NULLDOM_NOTDEF_MSG = """ This trait is not defined for NullDomain, ∅."""

expanse(::NullDomain) = throw(ArgumentError(NULLDOM_NOTDEF_MSG))
isperiodic(::NullDomain, args...) = throw(ArgumentError(NULLDOM_NOTDEF_MSG))
boundaries(::NullDomain) = ()
domain_tag(::NullDomain) = throw(ArgumentError(NULLDOM_NOTDEF_MSG))
boundary_tag(::NullDomain) = throw(ArgumentError(NULLDOM_NOTDEF_MSG))
bounding_box(::NullDomain) = throw(ArgumentError(NULLDOM_NOTDEF_MSG))

×(::NullDomain, ::AbstractDomain) = ∅
×(::AbstractDomain, ::NullDomain) = ∅

deform(::NullDomain, args...) = ∅

# Base.in(::Union{Number, AbstractDomain}, ::NullDomain) = false

###
# PointDomain
###

"""
Zero-dimensional domain representing a point.
"""
struct PointDomain{T<:Number} <: AbstractDomain{T, 0}
    x::T
    tag::Symbol

    function PointDomain(x::Number, tag::Union{Symbol, Nothing})
        tag = isnothing(tag) ? :NoTag : tag
        new{eltype(x)}(x, tag)
    end
end

PointDomain(x; tag = nothing) = PointDomain(x, tag)

function (::Type{T})(dom::PointDomain) where{T<:Number}
    @set! dom.x = T(dom.x)
end

POINTDOM_NOTDEF_MSG = """ This trait is not defined for PointDomain."""

expanse(::PointDomain{T}) where{T} = (zero(T),)
isperiodic(::PointDomain, args...) = throw(ArgumentError(POINTDOM_NOTDEF_MSG))
boundaries(::PointDomain) = (∅,)
domain_tag(dom::PointDomain) = dom.tag
function boundary_tag(::PointDomain, i)
    if i == 1 
        nothing
    else
        throw(ArgumentError("i > num_boundaries(::PointDomain)")) 
    end
end
bounding_box(dom::PointDomain) = dom

###
# IntervalDomain
###

"""
1D open interval containing `x` such that `x0 < x < x1`.
"""
struct IntervalDomain{T, Tp} <: AbstractDomain{T, 1}
    p0::Tp
    p1::Tp
    periodic::Bool
    tag::Symbol

    function IntervalDomain(p0::PointDomain, p1::PointDomain,
        periodic::Bool, tag::Union{Symbol, Nothing})

        @assert p0.x < p1.x "x0 < x1"
        T = promote_type(eltype.((p0, p1))...)
        p0 = T(p0)
        p1 = T(p1)
        tag = isnothing(tag) ? :NoTag : tag

        new{T, typeof(p0)}(p0, p1, periodic, tag)
    end
end

function IntervalDomain(x0, x1;
    periodic = false,
    tag = nothing,
    boundary_tags = (nothing, nothing),
)

    p0 = PointDomain(x0; tag = boundary_tags[1])
    p1 = PointDomain(x1; tag = boundary_tags[2])

    IntervalDomain(p0, p1, periodic, tag)
end

function (::Type{T})(int::IntervalDomain) where {T <: Number}
    IntervalDomain(T(int.p0), T(int.p1), int.periodic, int.tag)
end

expanse(dom::IntervalDomain) = (dom.p1.x - dom.p0.x,)
function isperiodic(dom::IntervalDomain, d::Integer)
    if d == 1
        dom.periodic
    else
        throw(ArgumentError("d > dims(dom)"))
    end
end
boundaries(dom::IntervalDomain) = (dom.p0, dom.p1,)
domain_tag(dom::IntervalDomain) = dom.tag
bounding_box(dom::IntervalDomain) = dom
#
