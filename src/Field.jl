#
# TODO: Fourier spectral, SpectralElements will require different functionality
#   - move broadcasting/indexing to AbstractField,
#   - if Fourier can be handled with vanilla Field, then rename it
#     SpectralField, else rename to LagrangePolyField
#
# TODO: create SpectralElementField and overload inner product by
#   - overload *(Adjoint{Field}, Field), norm(::Field, 2)
#
""" Scalar function field in D-dimensional space over a spectral basis"""
struct Field{T,D,Tarr <: AbstractArray{T,D}} <: AbstractField{T,D}
    array::Tarr
end

# display
function Base.summary(io::IO, u::Field{T,D}) where{T,D}
    println(io, "$(D)D scalar field of type $T")
    Base.show(io, typeof(u))
end
function Base.show(io::IO, ::MIME"text/plain", u::Field)
    iocontext = IOContext(io, :compact => true, :limit => true)
    Base.summary(iocontext, u)
    Base.show(iocontext, MIME"text/plain"(), u.array)
    println()
end

# vector indexing
Base.IndexStyle(::Field) = IndexLinear()
Base.getindex(u::Field, i::Int) = getindex(u.array, i)
Base.setindex!(u::Field, v, i::Int) = setindex!(u.array, v, i)
Base.size(u::Field) = (length(u.array),)

# allocation
function Base.similar(u::Field, ::Type{T}=eltype(u), dims::Dims=size(u.array)) where{T}
    Field(similar(u.array, T, dims))
end
Base.zero(u::Field, dims::Dims) = zero(u) # ignore dims since <: AbstractVector
function Field{T,D, Tarr}(::UndefInitializer, n) where{T,D,Tarr} # TODO - Krylov.jl hack. 
    vec = Vector{T}(undef, n)
    Field(vec)
end

# broadcast
Base.Broadcast.BroadcastStyle(::Type{<:Field}) = Broadcast.ArrayStyle{Field}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Field}},
                      ::Type{ElType}) where ElType
  u = find_fld(bc)
  Field(similar(Array{ElType}, axes(u.array)))
end

find_fld(bc::Base.Broadcast.Broadcasted) = find_fld(bc.args)
find_fld(args::Tuple) = find_fld(find_fld(args[1]), Base.tail(args))
find_fld(x) = x
find_fld(::Tuple{}) = nothing
find_fld(a::Field, rest) = a
find_fld(::Any, rest) = find_fld(rest)

LinearAlgebra.norm(u::Field, p::Real=2) = norm(u.array, p)

function LinearAlgebra.axpy!(a::Number, x::Field{<:Number,D}, y::Field{<:Number,D}) where{D}
    axpy!(a, x.array, y.array)
end

function LinearAlgebra.axpby!(a::Number, x::Field{<:Number,D}, b::Number, y::Field{<:Number,D}) where{D}
    axpby!(a, x.array, b, y.array)
end

function LinearAlgebra.dot(u::Field{<:Number,D}, v::Field{<:Number,D}) where{D}
    dot( _vec(u.array), _vec(v.array))
end

#TODO - implement lazy adjoint Adjoint(u::Field). overload ', adjoint. then
# overload (u' * A) \defeq  (A' * u)' in OperatorBasics.jl
#
#
