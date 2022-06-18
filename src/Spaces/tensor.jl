#
# do it for 2D first and nest later

struct TensorProductSpace{
                          T,
                          D1+D2,
                          I<:AbstractSpace{<:Number, D1},
                          O<:AbstractSpace{<:Number, D2},
                         } <: AbstractTensorProductSpace{T,D1+D2}
    inner::I
    outer::O
end

function gradOp(space::TensorProductSpace)
    @unpack inner, outer = space
    Di = gradOp(inner)
    Do = gradOp(outer)

    # make tensor product

#   [
#    Di
#    Do
#   ]
end

function massOp(space::TensorProductSpace)
    # multiply mass matrices
end
#
