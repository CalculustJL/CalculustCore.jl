#
###
# Product Domain
###

""" D-dimensional logically reectangular domain """
struct ProductDomain{T, D, Tdom} <: AbstractDomain{T, D}
    domains::Tdom

    function ProductDomain(domains)
        T = promote_type(eltype.(domains)...)
        D = sum(dims.(domains); init = 0)

        new{T, D, typeof(domains)}(domains)
    end
end

ProductDomain(domains::AbstractDomain...) = ProductDomain(domains)

function (::Type{T})(int::ProductDomain) where {T <: Number}
    ProductDomain(T.(int.domains)...)
end

×(doms::AbstractDomain...) = reduce(×, doms)
×(dom1::AbstractDomain, dom2::AbstractDomain) = ProductDomain(dom1, dom2)
×(dom1::ProductDomain, dom2::AbstractDomain) = ProductDomain(dom1.domains..., dom2)
×(dom1::AbstractDomain, dom2::ProductDomain) = ProductDomain(dom1, dom2.domains...)
×(dom1::ProductDomain, dom2::ProductDomain) = ProductDomain(dom1.domains..., dom2.domains...)

function expanse(dom::ProductDomain{T}) where{T}

    e = ()
    for d in dom.domains
        e = (e..., expanse(d)...)
    end

    T.(e)
end

function isperiodic(dom::ProductDomain{T,D}, dim::Integer) where{T,D}
    if dim ≤ D
        Ds = dims.(dom.domains)
        cs = cumsum(Ds)
        d = findfirst(n -> n ≥ dim, cs)

        isperiodic(dom.domains[d], dim + 1 - d)
    else
        throw(ArgumentError("d > dims(dom)"))
    end
end

function boundaries(dom::ProductDomain{T,D}) where{T,D}

    doms = dom.domains
    bdr = ()

    for d in eachindex(doms)
        _dom = doms[d]
        _bdr = boundaries(_dom)

        for i in eachindex(_bdr)
            # tag = domain_tag(_bdr[i])
            __bdr = ×(doms[begin:d-1]..., _bdr[i], doms[d+1:end]...)

        bdr = (bdr..., __bdr)
        end
    end

    bdr
end

function domain_tag(dom::ProductDomain)

    tags = domain_tag.(dom.domains)
    tags = unique(tags)

    all(isequal(:NoTag), tags) && return :NoTag

    # rm :NoTag
    i = findall(:NoTag, tags)
    tags = (tags[begin:i-1]..., tags[i+1:end])

    strs = string.(tag)
    str = strs[1]
    for i in 2:length(strs)
        str = string(str, " × ", strs[i])
    end

    Symbol(str)
end

bounding_box(dom::ProductDomain) = dom
#
