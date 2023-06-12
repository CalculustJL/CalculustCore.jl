#
"""
Tensor product of spaces. like extrusion
 do it for 2D first and nest later
"""
struct TensorProductSpace{
                          T,
                          D1 + D2,
                          I <: AbstractSpace{<:Any, D1},
                          O <: AbstractSpace{<:Any, D2}
                          } <: AbstractSpace{T, D1 + D2}
    inner::I
    outer::O
end

## TODO - add short circuits in downstream packages

#function domain(V::TensorProductSpace)
#    domain(outer) ⊗ domain(inner)
#end
#
#function points(V::TensorProductSpace)
#end
#
#function modes(V::TensorProductSpace)
#end
#
#function mass_matrix(V::TensorProductSpace)
#    DiagonalOperator(V.mass_matrix)
#end
#
#function global_numbering(V::TensorProductSpace)
#end
#
#function boundary_nodes()
#end

# vector calculus

function gradientOp(V::TensorProductSpace)
    @unpack inner, outer = V
    Di = gradientOp(inner)
    Do = gradientOp(outer)

    @error "method not implemented"

    # make tensor product

    #   [
    #    Di
    #    Do
    #   ]
end

function massOp(V::TensorProductSpace)
    # multiply mass matrices
    @error "method not implemented"
end
#
