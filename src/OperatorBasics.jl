#
# TODO
#   - caching doesn't work. fix it with mutable structs
#   - fuse composed operators if possible
#   - avoid caching for AffineOp() created with +,-,λ* where
#     the second op is just Identity, or Null
#
###
# AbstractOperator interface
###

SciMLBase.has_adjoint(::AbstractOperator) = true
SciMLBase.has_mul(::AbstractOperator) = true
SciMLBase.has_mul!(::AbstractOperator) = true

function (A::AbstractOperator{<:Number,D})(u::AbstractField{<:Number,D}, p, t::Number) where{D}
    A = SciMLBase.update_coefficients(A, u, p ,t)
    A * u
end

function (A::AbstractOperator{<:Number,D})(du::AbstractField{<:Number,D}, u, p, t::Number) where{D}
    SciMLBase.update_coefficients!(A, u, p ,t)
    mul!(du, A, u)
end

function LinearAlgebra.mul!(v::AbstractField{<:Number,D},A::AbstractOperator{<:Number,D},u::AbstractField{<:Number,D}) where{D}
    ArgumentError("LinearAlgebra.mul! not defined for $A")
end

function Base.:\(A::AbstractOperator{<:Number,D},u::AbstractField{<:Number,D}) where{D}
    ArgumentError("Operator inversion not defined for $A")
end

function Base.:*(A::AbstractOperator{<:Number,D}, B::AbstractOperator{<:Number,D}) where{D}
#   A ∘ B # https://github.com/vpuri3/PDEInterfaces.jl/issues/2
    ComposedOp(B, A)
end

# caching
function init_cache(A::AbstractOperator,u)
    @error "Caching behaviour not defined for $A"
end

function set_cache(A::AbstractOperator, cache)
    @set! A.cache = cache
    @set! A.isunset = false
    return A
end

Base.size(A::AbstractOperator, d::Integer) = d <= 2 ? size(A)[d] : 1

###
# NullOp
###

""" (Square) Null operator """
struct NullOp{D} <: AbstractOperator{Bool,D} end

Base.zero(::AbstractOperator{<:Number, D}) where{D} = NullOp{D}()
Base.zero(::Type{<:AbstractOperator{<:Number, D}}) where{D} = NullOp{D}()

Base.adjoint(Z::NullOp) = Z
issquare(::NullOp) = true

Base.:*(::NullOp{D}, u::AbstractField{<:Number,D}) where{D} = u * false

function LinearAlgebra.mul!(v::AbstractField{<:Number,D}, ::NullOp{D}, u::AbstractField{<:Number,D}) where{D}
    lmul!(false, v)
end

# overloads
for op in (
           :*, :∘,
          )
    @eval $op(::NullOp{D}, ::AbstractOperator{<:Number,D}) where{D} = NullOp{D}()
    @eval $op(::AbstractOperator{<:Number,D}, ::NullOp{D}) where{D} = NullOp{D}()
end

Base.:+(::NullOp{D}, A::AbstractOperator{<:Number,D}) where{D} = A
Base.:+(A::AbstractOperator{<:Number,D}, ::NullOp{D}) where{D} = A

Base.:-(::NullOp{D}, A::AbstractOperator{<:Number,D}) where{D} = -A
Base.:-(A::AbstractOperator{<:Number,D}, ::NullOp{D}) where{D} = A

###
# IdentityOp
###

""" (Square) Identity operator """
struct IdentityOp{D} <: AbstractOperator{Bool,D} end

Base.one(::AbstractOperator{<:Number, D}) where{D} = IdentityOp{D}()
Base.one(::Type{<:AbstractOperator{<:Number, D}}) where{D} = IdentityOp{D}()

Base.adjoint(Id::IdentityOp) = Id
Base.inv(Id::IdentityOp) = Id

issquare(::IdentityOp) = true

SciMLBase.has_ldiv(::IdentityOp) = true
SciMLBase.has_ldiv!(::IdentityOp) = true

Base.:*(::IdentityOp{D}, u::AbstractField{<:Number,D}) where{D} = copy(u)
Base.:\(::IdentityOp{D}, u::AbstractField{<:Number,D}) where{D} = copy(u)

function LinearAlgebra.mul!(v::AbstractField{<:Number,D}, ::IdentityOp{D}, u::AbstractField{<:Number,D}) where{D}
    copy!(v, u)
end

function LinearAlgebra.ldiv!(v::AbstractField{<:Number,D}, ::IdentityOp{D}, u::AbstractField{<:Number,D}) where{D}
    copy!(v, u)
end

function LinearAlgebra.ldiv!(Id::IdentityOp{D}, u::AbstractField{<:Number,D}) where{D}
    u
end

# overloads
for op in (
           :*, :∘,
          )
    @eval $op(::IdentityOp{D}, A::AbstractOperator{<:Number,D}) where{D} = A
    @eval $op(A::AbstractOperator{<:Number,D}, ::IdentityOp{D}) where{D} = A
end

Base.:/(A::AbstractOperator{<:Number,D}, ::IdentityOp{D}) where{D} = A

for op in (
           :*, :∘,
          )
    @eval $op(::NullOp{D}, ::IdentityOp{D}) where{D} = NullOp{D}()
    @eval $op(::IdentityOp{D}, ::NullOp{D}) where{D} = NullOp{D}()
end

###
# AffineOp
###

""" Lazy affine operator combinations αA + βB """
struct AffineOp{T,D,
                Ta <: AbstractOperator{<:Number,D},
                Tb <: AbstractOperator{<:Number,D},
                Tα,Tβ,Tc
               } <: AbstractOperator{T,D}
    A::Ta
    B::Tb
    α::Tα
    β::Tβ

    cache::Tc
    isunset::Bool

    function AffineOp(A::AbstractOperator{Ta,D}, B::AbstractOperator{Tb,D}, α, β,
                      cache = nothing, isunset = cache === nothing) where{Ta,Tb,D}
        T = promote_type(Ta,Tb)
        new{T,D,typeof(A),typeof(B),typeof(α),typeof(β),typeof(cache)}(A, B, α, β, cache, isunset)
    end
end

issquare(A::AffineOp) = issquare(A.A) & issquare(A.B)
function Base.adjoint(A::AffineOp)
    if issquare(A)
        AffineOp(A.A',A.B',A.α, A.β, A.cache, A.isunset)
    else
        AffineOp(A.A',A.B',A.α', A.β')
    end
end

function Base.convert(::Type{AffineOp}, A::AbstractOperator{<:Number,D}) where{D}
#   A
    AffineOp(A, NullOp{D}(), true, false)
end

function init_cache(A::AffineOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    cache = A.B * u
end

function Base.:*(A::AffineOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    @unpack A, B, α, β = A
    if iszero(α) | (A isa NullOp)
        β * (B * u)
    elseif iszero(β) | (B isa NullOp)
        α * (A * u)
    else
        α * (A * u) + β * (B * u)
    end
end

function LinearAlgebra.mul!(v::AbstractField{<:Number,D}, A::AffineOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    mul!(v, A.A, u)
    lmul!(A.α, v)

    if iszero(α) | (A isa NullOp)
        mul!(v, B, u)
        lmul!(β, v)
        return v
    elseif iszero(β) | (B isa NullOp)
        mul!(v, A, u)
        lmul!(α, v)
        return v
    end

    if isunset
        cache = init_cache(A,u)
        A = set_cache(A, cache)
    else
        mul!(A.cache, A.B, u)
    end

    lmul!(A.β, A.cache)
    axpy!(true, A.cache, v)
end

function Base.:+(A::AbstractOperator{<:Number,D}, B::AbstractOperator{<:Number,D}) where{D}
    AffineOp(A, B, true, true)
end

function Base.:-(A::AbstractOperator{<:Number,D}, B::AbstractOperator{<:Number,D}) where{D}
    AffineOp(A, B, true, -true)
end

function Base.:+(A::AbstractOperator{<:Number,D}, λ::Number) where{D}
    Id = IdentityOp{D}()
    AffineOp(A, Id, true, λ)
end

function Base.:+(λ::Number, A::AbstractOperator{<:Number,D}) where{D}
    Id = IdentityOp{D}()
    AffineOp(A, Id, true, λ)
end

function Base.:-(A::AbstractOperator{<:Number,D}, λ::Number) where{D}
    Id = IdentityOp{D}()
    AffineOp(A, Id, -true, λ)
end

function Base.:-(λ::Number, A::AbstractOperator{<:Number,D}) where{D}
    Id = IdentityOp{D}()
    AffineOp(Id, A, λ, -true)
end

function Base.:-(A::AbstractOperator{<:Number,D}) where{D}
    Z = NullOp{D}()
    AffineOp(A, -true, false, Z)
end

function Base.:*(A::AbstractOperator{<:Number,D}, λ::Number) where{D}
    Z = NullOp{D}()
    AffineOp(A, Z, λ, false)
end

function Base.:*(λ::Number, A::AbstractOperator{<:Number,D}) where{D}
    Z = NullOp{D}()
    AffineOp(A, Z, λ, false)
end

function Base.:/(A::AbstractOperator{<:Number,D}, λ::Number) where{D}
    Z = NullOp{D}()
    AffineOp(A, Z, -true, λ)
end

###
# ComposedOp
###

""" Lazy Composition A ∘ B """
struct ComposedOp{T,D,Ti,To,Tc} <: AbstractOperator{T,D}
    inner::Ti
    outer::To

    cache::Tc
    isunset::Bool

    function ComposedOp(inner::AbstractOperator{Ti,D},
                        outer::AbstractOperator{To,D},
                        cache = nothing,
                        isunset::Bool = cache === nothing
                       ) where{Ti,To,D}
#       @assert size(outer, 1) == size(inner, 2)
        T = promote_type(Ti, To)
        isunset = cache === nothing
        new{T,D,typeof(inner),typeof(outer),typeof(cache)}(inner, outer, cache, isunset)
    end
end

function Base.:∘(outer::AbstractOperator{Number,D}, inner::AbstractOperator{Number,D}) where{D}
    ComposedOp(inner, outer)
end

Base.size(A::ComposedOp) = (size(A.outer, 1), size(A.inner, 2))
# https://github.com/vpuri3/PDEInterfaces.jl/issues/2
Base.adjoint(A::ComposedOp) = ComposedOp(A.outer', A.inner') # A.inner' ∘ A.outer'
Base.inv(A::ComposedOp) = ComposedOp(inv(A.outer), inv(A.inner)) # inv(A.inner) ∘ inv(A.outer)

SciMLBase.has_ldiv(A::ComposedOp) = has_ldiv(A.inner) & has_ldiv(A.outer)
SciMLBase.has_ldiv!(A::ComposedOp) = has_ldiv!(A.inner) & has_ldiv!(A.outer)
issquare(A::ComposedOp) = issquare(A.inner) & issquare(A.outer)

function init_cache(A::ComposedOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    cache = A.inner * u
    return cache
end

function Base.:*(A::ComposedOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    @unpack inner, outer = A
    outer * (inner * u)
end

function Base.:\(A::ComposedOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    @assert has_ldiv(A)
    @unpack inner, outer = A

    outer \ (inner \ u)
end

function LinearAlgebra.mul!(v::AbstractField{<:Number,D}, A::ComposedOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    @unpack inner, outer = A

    if A.isunset
        cache = init_cache(A, u)
        A = set_cache(A, cache) # pass this A to calling context
    else
        mul!(cache, inner, u)
    end

    mul!(v, outer, cache)
end

function LinearAlgebra.ldiv!(A::ComposedOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    @assert has_ldiv!(A)
    @unpack inner, outer = A

    ldiv!(inner, u)
    ldiv!(outer, u)
end

function LinearAlgebra.ldiv!(v::AbstractField{<:Number,D}, A::ComposedOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    @assert has_ldiv!(A)
    @unpack inner, outer = A

    ldiv!(v, inner, u)
    ldiv!(outer, v)
end

###
# InverseOp
###

""" Lazy Inverse Operator """
struct InverseOp{T,D,Ta} <: AbstractOperator{T,D}
    A::Ta

    function InverseOp(A::AbstractOperator{T,D}) where{T,D}
        @assert issquare(A)
        new{T,D,typeof(A)}(A)
    end
end

function Base.:*(A::InverseOp{<:Number,D}, u::AbstractField{<:Number,D}) where{D}
    A.A \ u
end

SciMLBase.has_ldiv(::InverseOp) = true
SciMLBase.has_ldiv!(::InverseOp) = true
issquare(::InverseOp) = true

Base.inv(A::AbstractOperator) = InverseOp(A)
Base.size(A::InverseOp) = size(A.A)
Base.adjoint(A::InverseOp) = inv(A.A')
LinearAlgebra.ldiv!(y, A::InverseOp, x) = mul!(y, A.A, x)
LinearAlgebra.mul!(y, A::InverseOp, x) = ldiv!(y, A.A, x)

# overload inverse for matrices of operators
#function Base.inv(A::Matrix{<:AbstractOperator})
#end
#
