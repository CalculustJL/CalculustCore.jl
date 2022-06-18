#
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
        x1 = a1 + λ1 × x1(r1), ...,
        xD = aD + λD × xD(rD)
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


#""" D-Dimensional domain maps """
abstract type AbstractMap{D} end

# TODO - make a struct for mappings
struct Deformations{D,Tmap} <: AbstractMap{D}
  mapping::Tmap
  isrescaling::Bool
  isseparable::Bool
  isinvertibe::Bool # default to true
end

# AffineMap y = Ax + b
# LinearMap y = Ax
# SeparableMap: y = D x + b
# 

###
# Transfinite Interpolation
###

function gordon_hall end #TODO

"""
 Gordon Hall map - Transfinite interpolation
"""
function mapBoxes(box1::BoxDomain{<:Number,D}, box2::BoxDomain{<:Number,D}) where{D}
    mapping = nothing
end

function map_from_ref(domain, ref_domain;D=D) # TODO

    if domains_match(domain, ref_domain)
        return domain
    end

    xe = domain.endpoints

    mapping = domain.mapping ==! nothing ? domain.mapping : identity#(r -> r)
    mapping = mapping ∘ map
end

#=
function fill_box_with_space(dom1::BoxDomain{<:Number,D},
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

function polarMap(r, θ)
    x = @. r * cos(θ)
    y = @. r * sin(θ)
    x, y
end

function reference_box(D; periodic_dirs=())
    domain = BoxDomain()
    for i=1:D
        periodic = i ∈ periodic_dirs
        tag1 = Symbol("Lower$(i)")
        tag2 = Symbol("Upper$(i)")
        tags = (tag1, tag2)
        interval = IntervalDomain(-true, true; periodic=periodic, bdry_tags=tags)
        domain *= interval
    end
    domain
end

function annulus_2D(r0, r1)
    intR = IntervalDomain(r0, r1, false, (:Inner, :Outer))
    intθ = IntervalDomain(-π,  π, true , (:Periodic, :Periodic))

    dom = intR * intθ

    deform(dom, polar)
end
