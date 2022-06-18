#
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

