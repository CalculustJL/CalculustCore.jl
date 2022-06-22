
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
Base.size(space::TransformedSpace) = size(space.space)

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

###
# vector calculus - modal space <-> modal space
###

function gradOp(space::TransformedSpace) # ∇
end

#
