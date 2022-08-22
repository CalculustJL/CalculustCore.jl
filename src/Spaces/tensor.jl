#
"""
Tensor product of spaces. like extrusion
 do it for 2D first and nest later
"""
struct TensorProductSpace{
                          T,
                          D1+D2,
                          I<:AbstractSpace{<:Any, D1},
                          O<:AbstractSpace{<:Any, D2},
                         } <: AbstractSpace{T,D1+D2}
    inner::I
    outer::O
end

function domain(space::TensorProductSpace)
    domain(outer) ⊗ domain(inner)
end

function points(space::TensorProductSpace)
end

function quadratures(space::TensorProductSpace)
    (quadratures(inner)..., quadratures(outer)...)
end

function mass_matrix(space::TensorProductSpace)
    DiagonalOperator(space.mass_matrix)
end

function local_numbering(space::TensorProductSpace)
end

function boundary_nodes()
end

# vector calculus

function gradientOp(space::TensorProductSpace)
    @unpack inner, outer = space
    Di = gradientOp(inner)
    Do = gradientOp(outer)

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
