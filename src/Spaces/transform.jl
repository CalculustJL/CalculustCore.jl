#
"""
Lazy wrapper for transforming

    physical space -> modal space
"""
struct TransformedSpace{T, D, Tsp <: AbstractSpace{T, D}} <: AbstractSpace{T, D}
    space::Tsp
end

transform(V::AbstractSpace) = TransformedSpace(V)
transform(Vh::TransformedSpace) = Vh.space

###
# interface
###
Base.size(Vh::TransformedSpace) = mode_size(Vh.space)

domain(Vh::TransformedSpace) = domain(Vh.space)
points(Vh::TransformedSpace) = modes(Vh.space)
modes(Vh::TransformedSpace) = points(Vh.space)

function Domains.deform(Vh::TransformedSpace, args...)

    V = transform(Vh)
    V = deform(V, args...)

    transform(V)
end

function form_transform(Vh::TransformedSpace, u::AbstractArray; kwargs...)
    F = transformOp(Vh)
    V = form_transform(Vh, F * u; kwargs...)

    transform(V)
end

#function Domains.⊗(Vh::TransformedSpace, W::AbstractSpace)
#end

###
# vector calculus - modal space <-> modal space
###

#function massOp(Vh::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function gradientOp(Vh::TransformedSpace)
#end
#
#function hessianOp(Vh::TransformedSpace)
#end
#
#function laplaceOp(Vh::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function diffusionOp(ν::Number, Vh::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function diffusionOp(ν::AbstractVector, Vh::TransformedSpace, discr::AbstractDiscretization)
#end
#
#function advectionOp(vel::NTuple{D}, Vh::TransformedSpace{<:Any,D}) where{D}
#end

function transformOp(Vh::TransformedSpace)
    ftr = transformOp(Vh.space)
    inv(ftr)
end
#
