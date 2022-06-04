#
###
# Matrix Operator
###

""" 1D mat-vecs """
struct MatrixOp{T,Tm<:AbstractMatrix{T}} <: AbstractOperator{T,1}
    mat::Tm
end

@forward MatrixOp.mat (
                       issquare, SciMLBase.has_ldiv, SciMLBase.has_ldiv!
                      )

Base.size(A::MatrixOp) = size(A.mat)
Base.adjoint(A::MatrixOp) =  MatrixOp(A.mat')
Base.inv(A::MatrixOp) = MatrixOp(inv(A.mat))

function Base.:*(A::MatrixOp, u::AbstractField{<:Number,1})
    A.mat * _vec(u)
end

function Base.:\(A::MatrixOp, u::AbstractField{<:Number,1})
    A.mat \ _vec(u)
end

function LinearAlgebra.mul!(v::AbstractField{<:Number,1}, A::MatrixOp, u::AbstractField{<:Number,1})
    mul!(_vec(v), A.mat, _vec(u))
    return v
end

function LinearAlgebra.ldiv!(v::AbstractField{<:Number,1}, A::MatrixOp, u::AbstractField{<:Number,1})
    ldiv!(_vec(v), A.mat, _vec(u))
    return v
end

# fusion
for op in (
           :+, :-, :*,
          )
    @eval function Base.$op(A::MatrixOp, B::MatrixOp)
        M = A.mat * B.mat
        MatrixOp(M)
    end
end

###
# Diagonal Operator
###

""" Diagonal Scaling Operator """
struct DiagonalOp{T,D,Tdiag<:AbstractField{T,D}} <: AbstractOperator{T,D}
    diag::Tdiag # diagonal vector
end

Base.size(A::DiagonalOp) = size(Diagonal(A.diag))
Base.adjoint(A::DiagonalOp) = A
Base.inv(A::DiagonalOp) = DiagonalOp(1 ./ A.diag)

SciMLBase.has_ldiv(::DiagonalOp) = true
SciMLBase.has_ldiv!(::DiagonalOp) = true
issquare(::DiagonalOp) = true

function Base.:*(A::DiagonalOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    Diagonal(A.diag) * _vec(u)
end

function Base.:\(A::DiagonalOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    Diagonal(A.diag) \ _vec(u)
end

function LinearAlgebra.mul!(v::AbstractField{<:Number,D}, A::DiagonalOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    mul!(_vec(v), Diagonal(A.diag), _vec(u))
    return v
end

function LinearAlgebra.ldiv!(v::AbstractField{<:Number,D}, A::DiagonalOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    ldiv!(_vec(v), Diagonal(A.diag), _vec(u))
    return v
end

function LinearAlgebra.ldiv!(A::DiagonalOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    ldiv!(Diagonal(A.diag), _vec(u))
    return u
end

# fusion
for op in (
           :+, :-, :*,
          )
    @eval function Base.$op(A::DiagonalOp, B::DiagonalOp)
        Diag = $op(Diagonal(A.diag), Diagonal(B.diag))
        DiagonalOp(Diag.diag)
    end

    @eval function Base.$op(λ::Number, A::DiagonalOp)
        diag = $op(λ, A.diag)
        DiagonalOp(diag)
    end

    @eval function Base.$op(A::DiagonalOp, λ::Number)
        diag = $op(A.diag, λ)
        DiagonalOp(diag)
    end
end

function Base.:*(A::MatrixOp, D::DiagonalOp{<:Number,1})
    M = A.mat * Diagonal(D.diag)
    MatrixOp(M)
end

function Base.:*(D::DiagonalOp{<:Number,1}, A::MatrixOp)
    M = Diagonal(D.diag) * A.mat
    MatrixOp(M)
end

function Base.:/(A::DiagonalOp, λ::Number)
    diag = A.diag ./ λ
    DiagonalOp(diag)
end

function Base.:/(A::DiagonalOp, B::DiagonalOp)
    Diag = Diagonal(A.diag) / Diagonal(B.diag)
    DiagonalOp(Diag.diag)
end

function Base.:\(A::DiagonalOp, B::DiagonalOp)
    Diag = Diagonal(A.diag) \ Diagonal(B.diag)
    DiagonalOp(Diag.diag)
end

LinearAlgebra.:lmul!(a::Number,B::DiagonalOp) = lmul!(a,B.diag)
LinearAlgebra.:rmul!(A::DiagonalOp,b::Number) = rmul!(A.diag,b)

###
# Tensor Product Operators
###

"""
Tensor product operator

(B ⊗ A) * u
"""
function tensor_product!(V,U,A,B,cache) # 2D
    """ V .= A * U * B' """
    mul!(cache, A, U)
    mul!(V, cache, B')
end

"""
Tensor product operator

(C ⊗ B ⊗ A) * u
"""
function tensor_product!(V,U,A,B,Ct,cache1,cache2) # 3D
    szU = size(U)
    U_re = _reshape(U, (szU[1], szU[2]*szU[3]))
    mul!(cache1, A, U_re)

    # TODO
    # B op - write to cache2. use views

    szC = size(cache2)
    C_re = _reshape(cache2, (szC[1]*szC[2], szC[3]))
    mul!(V, C_re, Ct')

    V
end

###
# TensorProductOp2D
###

""" 2D Tensor Product Operator """
struct TensorProductOp2D{T,Ta,Tb,Tc} <: AbstractTensorProductOperator{T,2}
    A::Ta
    B::Tb

    cache::Tc
    isunset::Bool

    function TensorProductOp2D(A, B, cache = nothing, isunset = cache === nothing)
        T = promote_type(eltype(A), eltype(B))
        new{T,typeof(A),typeof(B),typeof(cache)}(A, B, cache, isunset)
    end
end

Base.size(A::TensorProductOp2D) = size(A.A) .* size(A.B)
issquare(A::TensorProductOp2D) = issquare(A.A) & issquare(A.B)

function Base.adjoint(A::TensorProductOp2D)
    if issquare(A)
        TensorProductOp2D(A.A', A.B', A.cache)
    else
        TensorProductOp2D(A.A', A.B')
    end
end

function Base.:*(A::TensorProductOp2D{<:Number}, u::AbstractField{<:Number,D}) where{D}
    v = copy(u)
    @set! v.array = A.A * u.array * A.B' # get type information from u
end

function init_cache(A::TensorProductOp2D, U)
    cache = A.A * U
end

function LinearAlgebra.mul!(v::AbstractField{<:Number,2}, A::TensorProductOp2D, u::AbstractField{<:Number,2})
    U = u.array
    V = v.array

    if A.isunset
        cache = init_cache(A, U)
        A = set_cache(A, cache)
        mul!(V, cache, A.B')
        return v
    end
    
    tensor_product!(V, U, A.A, A.B, A.cache)
    return v
end

# fusion
for op in (
           :+ , :- , :* , :/, :\,
          )
    @eval function Base.$op(A::TensorProductOp2D, B::TensorProductOp2D)
        A = $op(A.A, B.A)
        B = $op(A.B, B.B)
        TensorProductOp2D(A, B)
    end
end
#
