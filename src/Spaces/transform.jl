
# TODO allow for evolution in either transformed space, or physical space
struct TransformedSpace
    space
end

function transform(space::AbstractSpace)
    TransformedSpace(space)
end

function transform(space::TransformedSpace)
    space.space
end
#
