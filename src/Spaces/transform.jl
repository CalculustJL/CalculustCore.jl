
"""
 physical space -> modal space
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
points(space::TransformedSpace) = points(space.space)
modes(space::TransformedSpace) = modes(space.space)
quadratures(space::TransformedSpace) = quadratures(space.space)

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

function massOp(space::TransformedSpace, discr::AbstractDiscretization) # ∫
end

function gradOp(space::TransformedSpace) # ∇
end

function hessianOp(space::TransformedSpace) # ∇²
end

function laplaceOp(space::TransformedSpace, discr::AbstractDiscretization) # Δ
end

function diffusionOp(ν::Number, space::TransformedSpace, discr::AbstractDiscretization) # νΔ
end

function diffusionOp(ν::AbstractVector, space::TransformedSpace, discr::AbstractDiscretization) # ∇⋅(ν∇)
end

function advectionOp(vel::NTuple{D}, space::TransformedSpace{<:Any,D}) where{D}
end

#
