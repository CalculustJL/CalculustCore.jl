#
""" use Base.ReshapedArray """
_reshape(a,dims::NTuple{D,Int}) where{D} = reshape(a,dims)
_reshape(a::Array, dims::NTuple{D,Int}) where{D} = Base.ReshapedArray(a, dims, ())

_vec(a) = vec(a)
_vec(a::AbstractVector) = a
_vec(a::AbstractArray) = _reshape(a,(length(a),))

""" check if operator(s) is square """
issquare(::UniformScaling) = true
issquare(A) = size(A,1) === size(A,2)
issquare(A...) = @. (&)(issquare(A)...)
#

# make this multidimensional by using the multidimensional indexing
# trick in domains
struct TensorProduct2DOperator{T} <: SciMLOperators.AbstractSciMLOperator{T}
    A
    B

    cache
    isset

    function TensorProduct2DOperator(A, B, cache, isset)
        T = promote_type(eltype.((A, B))...)
        isset = cache !== nothing
        new{T}(A, B, cache, isset)
    end
end

function TensorProduct2DOperator(A::AbstractMatrix, B::AbstractMatrix; cache = nothing)
    isset = cache !== nothing
    TensorProduct2DOperator(A, B, cache, isset)
end

Base.size(L::TensorProduct2DOperator) = size(A.A) .* size(A.B)
function Base.adjoint(L::TensorProduct2DOperator)
    TensorProduct2DOperator(L.A', L.B'; cache = issquare(A) ? L.cache : nothing)
end

function Base.:*(L::TensorProduct2DOperator, u::AbstractVector)
    sz = (size(L.A, 2), size(L.B, 2))
    u = _reshape(u, sz)
    v = L.A * u * L.B'
end
