#
###
# Product Domain
###

"""
$TYPEDEF
"""
struct ProductDomain{T, D, Tdom} <: AbstractDomain{T, D}
    domains::Tdom
    tag::Symbol

    function ProductDomain(domains, tag::Union{Symbol, Nothing})
        tag = isnothing(tag) ? :NoTag : tag
        T = promote_type(eltype.(domains)...)
        D = sum(ndims.(domains); init = 0)

        new{T, D, typeof(domains)}(domains, tag)
    end
end

function ProductDomain(domains::AbstractDomain...; tag = nothing)
    ProductDomain(domains, tag)
end

function (::Type{T})(int::ProductDomain) where {T <: Number}
    ProductDomain(T.(int.domains)...)
end

"""
Take cartesian product of `AbstractDomain`s `doms`

$SIGNATURES
"""
×(doms::AbstractDomain...) = reduce(×, doms)
×(dom1::AbstractDomain, dom2::AbstractDomain) = ProductDomain(dom1, dom2)
×(dom1::ProductDomain, dom2::AbstractDomain) = ProductDomain(dom1.domains..., dom2)
×(dom1::AbstractDomain, dom2::ProductDomain) = ProductDomain(dom1, dom2.domains...)
×(dom1::ProductDomain, dom2::ProductDomain) = ProductDomain(dom1.domains..., dom2.domains...)

function bounding_box(dom::ProductDomain)
    ProductDomain(bounding_box.(dom.domains)...)
end

function expanse(dom::ProductDomain{T}) where{T}
    e = ()
    for d in dom.domains
        e = (e..., expanse(d)...)
    end

    T.(e)
end

function isperiodic(dom::ProductDomain{T,D}) where{T,D}
    p = ()
    for d in dom.domains
        p = (p..., isperiodic(d)...)
    end

    p
end

function boundaries(dom::ProductDomain{T,D}) where{T,D}

    doms = dom.domains
    bdr = ()

    for d in eachindex(doms)
        _dom = doms[d]
        _bdr = boundaries(_dom)

        for i in eachindex(_bdr)
            __bdr = ×(doms[begin:d-1]..., _bdr[i], doms[d+1:end]...)

            bdr = (bdr..., __bdr)
        end
    end

    bdr
end

function domain_tag(dom::ProductDomain)

    if !tag_notdef(dom.tag)
        return dom.tag
    end

    tags = domain_tag.(dom.domains)
    tags = unique(tags)

    all(tag_notdef, tags) && return :NoTag

    # rm :NoTag
    i = findall(tag_notdef, tags)
    if length(i) === 1
        i = first(i)
        tags = (tags[begin:i-1]..., tags[i+1:end]...)
    end

    length(tags) == 1 && return first(tags)

    strs = string.(tags)
    str = strs[1]

    for i in 2:length(strs)
        str = string(str, " × ", strs[i])
    end

    Symbol(str)
end

function Base.show(io::IO, dom::ProductDomain)
    doms = dom.domains

    if length(doms) == 0
        print(io, "ProductDomain()")
        return
    end

    show(io, doms[1])
    for d in 2:length(doms)
        print(io, " × ")
        show(io, doms[d])
    end
end
#
