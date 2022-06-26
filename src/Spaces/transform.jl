#
"""
Lazy wrapper for transforming

    physical space -> modal space

Returned operators act on vectors in modal space. All other behaviour is the same.
"""
struct TransformedSpace{T,D,Tsp<:AbstractSpace{T,D}} <: AbstractSpace{T,D}
    space::Tsp
end

function transform(space::AbstractSpace)
    TransformedSpace(space)
end

function transform(space::TransformedSpace)
    space.space
end

###
# interface
###
Base.size(space::TransformedSpace) = size(modes(space))

domain(space::TransformedSpace) = domain(space.space)
points(space::TransformedSpace) = modes(space.space)
modes(space::TransformedSpace) = points(space.space)

function deform(space::TransformedSpace)
    phys = transform(space)
    def  = deform(space)

    transform(def)
end

#function Domains.⊗(space::TransformedSpace, space::AbstractSpace)
#end

###
# vector calculus - modal space <-> modal space
###

function massOp(space::TransformedSpace, discr::AbstractDiscretization)
end

function gradientOp(space::TransformedSpace)
end

function hessianOp(space::TransformedSpace)
end

function laplaceOp(space::TransformedSpace, discr::AbstractDiscretization)
end

function diffusionOp(ν::Number, space::TransformedSpace, discr::AbstractDiscretization)
end

function diffusionOp(ν::AbstractVector, space::TransformedSpace, discr::AbstractDiscretization)
end

function advectionOp(vel::NTuple{D}, space::TransformedSpace{<:Any,D}) where{D}
end

function transformOp(space::TransformedSpace)
    ftr = transformOp(space.space)
    inv(ftr)
end
#
