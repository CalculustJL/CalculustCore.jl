
# TODO allow for evolution in either transformed space, or physical space
struct TransformedSpace{T,D,Tsp<:AbstractSpace{T,D}} <: AbstractSpace{T,D}
    space::Tsp
end

function transform(space::AbstractSpace)
    TransformedSpace(space)
end

function transform(space::TransformedSpace)
    space.space
end

function domain(space::TransformedSpace)
end

function gradOp(space::TransformedSpace)
end
#
