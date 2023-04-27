#
"""
Lazy wrapper for transforming

    physical space -> modal space
"""
struct TransformedSpace{T, D, Tsp <: AbstractSpace{T, D}} <: AbstractSpace{T, D}
    space::Tsp
end

transform(space::AbstractSpace) = TransformedSpace(space)
transform(space::TransformedSpace) = space.space

###
# interface
###
Base.size(space::TransformedSpace) = mode_size(space.space)

domain(space::TransformedSpace) = domain(space.space)
points(space::TransformedSpace) = modes(space.space)
modes(space::TransformedSpace) = points(space.space)

function Domains.deform(space::TransformedSpace)
    phys = transform(space)
    def = deform(space)

    transform(def)
end

function form_transform(space::TransformedSpace, u::AbstractArray; kwargs...)
    F = transformOp(space)
    phys = form_transform(space, F * û; kwargs...)

    transform(phys)
end

#function Domains.⊗(space::TransformedSpace, space::AbstractSpace)
#end

###
# vector calculus - modal space <-> modal space
###

#function massOp(space::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function gradientOp(space::TransformedSpace)
#end
#
#function hessianOp(space::TransformedSpace)
#end
#
#function laplaceOp(space::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function diffusionOp(ν::Number, space::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function diffusionOp(ν::AbstractVector, space::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function advectionOp(vel::NTuple{D}, space::TransformedSpace{<:Any,D}) where{D}
#end

function transformOp(space::TransformedSpace)
    ftr = transformOp(space.space)
    inv(ftr)
end
#
