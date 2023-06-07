#
###
# DomainMap
###
"""
# Arguments
- reference domain

- mapping function
    (x1,...,xD) = map(r1, ..., rD)

- isrescaling - do optimizations if mapping is a simple rescaling
    x1 = a1 + λ1 × x1(r1),
    ...,
    xD = aD + λD × xD(rD)

-isseparable - do optimizations if mapping is separable
    x1 = x1(r1),
    ...,
    xD = xD(rD)

TODO - make types for AffineMap, LinearMap, SeparableMap, Transation, Rotation
TODO - write constructor that checks for method
"""
struct DomainMap{Tm}
    map::Tm
    isseparable::Bool
    isrescaling::Bool
end

function DomainMap(map; isseparable = false, isrescaling = false)
    DomainMap(map, isseparable, isrescaling)
end

###
# MappedDomain
###

"""
Deform D-dimensional domain via mapping
"""
struct MappedDomain{T, D, TD, TM} <: AbstractDomain{T, D} where{TD <: AbstractDomain{T, D}, TM <: DomainMap}
    domain::TD
    map::TM
end

function deform(dom::AbstractDomain, map::DomainMap)
    MappedDomain(dom, map)
end

identity_map(xyz::AbstractArray...) = xyz

function deform(dom::AbstractDomain, map = nothing;
                isseparable = false, isrescaling = false)

    map = if isnothing(map)
        IdentityMap()
    else
        _map = DomainMap(map; isseparable = isseparable,
                         isrescaling = isrescaling)
    end

    deform(dom, _map)
end

function (::Type{T})(dom::MappedDomain) where {T <: Number}
    MappedDomain(T(dom.domain), dom.map)
end

###
# Transfinite Interpolation
###

# function gordon_hall end #TODO

# """
#  Gordon Hall map - Transfinite interpolation
# """
# function mapBoxes(box1::BoxDomain{<:Number, D}, box2::BoxDomain{<:Number, D}) where {D}
#     mapping = nothing
# end
#
#=
# TODO Gordon Hall
# function mesh_domain(dom1::BoxDomain{<:Number,D},
#                      space::AbstractSpace{<:Number,D}
#                     ) where{D}
#
#     #function gordonHall(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zr,zs)
#     ze  = [-1,1]
#     Jer = interpMat(zr,ze)
#     Jes = interpMat(zs,ze)
#
#     xv = [xrm[1]   xrp[1]
#           xrm[end] xrp[end]]
#
#     yv = [yrm[1]   yrp[1]
#           yrm[end] yrp[end]]
#
#     xv = ABu(Jes,Jer,xv)
#     yv = ABu(Jes,Jer,yv)
#
#     #display(mesh(xv,yv,0*xv,0,90))
#
#     x = ABu([],Jer,vcat(xrm',xrp')) .+ ABu(Jes,[],hcat(xsm,xsp)) .- xv
#     y = ABu([],Jer,vcat(yrm',yrp')) .+ ABu(Jes,[],hcat(ysm,ysp)) .- yv
#
#     return
# end
=#
