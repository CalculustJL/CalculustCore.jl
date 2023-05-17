#=
# TODO Gordon Hall
function mesh_domain(dom1::BoxDomain{<:Number,D},
                     space::AbstractSpace{<:Number,D}
                    ) where{D}

    #function gordonHall(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zr,zs)
    ze  = [-1,1]
    Jer = interpMat(zr,ze)
    Jes = interpMat(zs,ze)

    xv = [xrm[1]   xrp[1]
          xrm[end] xrp[end]]

    yv = [yrm[1]   yrp[1]
          yrm[end] yrp[end]]

    xv = ABu(Jes,Jer,xv)
    yv = ABu(Jes,Jer,yv)

    #display(mesh(xv,yv,0*xv,0,90))

    x = ABu([],Jer,vcat(xrm',xrp')) .+ ABu(Jes,[],hcat(xsm,xsp)) .- xv
    y = ABu([],Jer,vcat(yrm',yrp')) .+ ABu(Jes,[],hcat(ysm,ysp)) .- yv

    return
end
=#

###
# Conveniences
###

function periodic_interval_tags(dim::Integer)
    tag = Symbol("D$(dim)_periodic")

    (tag, tag)
end

function default_interval_tags(dim::Integer)
    tag1 = Symbol("D$(dim)_inf")
    tag2 = Symbol("D$(dim)_sup")

    (tag1, tag2)
end

function RectangleDomain(x0, x1, y0, y1)
end

UnitIntervalDomain(; kwargs...) = UnitDomain(1; kwargs...)
UnitSquareDomain(; kwargs...) = UnitDomain(2; kwargs...)
UnitCubeDomain(; kwargs...) = UnitDomain(3; kwargs...)

UnitDomain(D; kwargs...) = HyperSquareDomain(D, 0.0, 1.0; kwargs...)
ChebyshevDomain(D; kwargs...) = HyperSquareDomain(D, -1.0, 1.0; kwargs...)
FourierDomain(D; kwargs...) = HyperSquareDomain(D, -π, π; periodic_dirs = 1:D, kwargs...)

function HyperSquareDomain(D::Integer, x0::Number, x1::Number; periodic_dirs = ())
    domain = ProductDomain()
    endpts = (x0, x1)

    for d in 1:D
        periodic = d ∈ periodic_dirs
        bdr_tags = periodic ? periodic_interval_tags(d) : default_interval_tags(d)

        interval = IntervalDomain(endpts...; periodic = periodic,
            tag = :Interior, boundary_tags = bdr_tags)

        domain = domain × interval
    end

    domain
end

function HyperRectangleDomain(D::Integer, endpoints::NTuple{2,<:Number}...;
    periodic_dirs = ())
    domain = ProductDomain()

    @assert length(endpoints) == D

    for d in 1:D
        periodic = d ∈ periodic_dirs
        endpts = endpoints[d]
        bdr_tags = periodic ? periodic_interval_tags(d) : default_interval_tags(d)

        interval = IntervalDomain(endpts...; periodic = periodic,
            tag = :Interior, boundary_tags = bdr_tags)

        domain = domain × interval
    end
end

function PolarMap()
    function polar(r, θ)
        x = r * cos(θ)
        y = r * sin(θ)
        x, y
    end

    DomainMap(polar)
end

function AnnulusDomain(rIn, rOut)
    intR = IntervalDomain(rIn, rOut, false, (:Inner, :Outer))
    intθ = IntervalDomain(-π, π, true, (:Periodic, :Periodic))

    dom = intR × intθ

    deform(dom, PolarMap())
end
#
