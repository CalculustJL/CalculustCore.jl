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

UnitBoxDomain(D; kwargs...) = FixedEndpointBoxDomain(D, 0.0, 1.0; kwargs...)
ChebyshevDomain(D; kwargs...) = FixedEndpointBoxDomain(D, -1.0, 1.0; kwargs...)
FourierDomain(D; kwargs...) = FixedEndpointBoxDomain(D, -π, π; periodic_dims = 1:D, kwargs...)

UnitIntervalDomain(; kwargs...) = UnitBoxDomain(1; kwargs...)
UnitSquareDomain(; kwargs...) = UnitBoxDomain(2; kwargs...)
UnitCubeDomain(; kwargs...) = UnitBoxDomain(3; kwargs...)

"""
`D`-dimensional `BoxDomain` with same endpoints for all dimensions.
"""
function FixedEndpointBoxDomain(D::Integer, x0::Number, x1::Number; kwargs...)
    endpoints = ()

    for _ in 1:D
        endpoints = (endpoints..., x0, x1)
    end

    BoxDomain(endpoints...; kwargs...)
end

"""
Generate a `D`-dimensional box domain with endpoints given by `endpoints`.
"""
function BoxDomain(endpoints::Real...;
                   periodic_dims = (),
                   tag = nothing,
                   boundary_tags = nothing)

    @assert length(endpoints) |> iseven "Number of endpoints must be even!"
    D = div(length(endpoints), 2)
    tag = isnothing(tag) ? :Interior : tag

    intervals = ()
    for d in 1:D
        idx = [2d - 1, 2d]

        periodic = d ∈ periodic_dims
        bdr_tags = if isnothing(boundary_tags)
            periodic ? periodic_interval_tags(d) : default_interval_tags(d)
        else
            boundary_tags[idx]
        end

        interval = IntervalDomain(endpoints[idx]...;
                                  periodic = periodic,
                                  boundary_tags = bdr_tags)

        intervals = (intervals..., interval)
    end

    ProductDomain(intervals...; tag = tag)
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
