#
Base.eltype(::AbstractSpace{T}) where{T} = T

"""
Dimension of underlying domain
"""
Domains.dims(::AbstractSpace{<:Any,D}) where{D} = D
Domains.dims(::AbstractArray{<:Any,D}) where{D} = D

function SciMLOperators.IdentityOperator(space::AbstractSpace)
    N = length(space)
    SciMLOperators.IdentityOperator{N}()
end

function SciMLOperators.NullOperator(space::AbstractSpace)
    N = length(space)
    SciMLOperators.NullOperator{N}()
end

"""
get number of points
"""
Base.size

"""
length of vector in space
"""
Base.length(space::AbstractSpace) = prod(size(space))

function Base.summary(io::IO, space::AbstractSpace{T,D}) where{T,D}
    type_color, no_color = SciMLBase.get_colorizers(io)
    print(io,
          type_color, nameof(typeof(space)),
          no_color," over domain ",
          type_color,typeof(domain(space)),
          no_color," with uType ",
          type_color,typeof(first(grid(space))),
          no_color
         )
end

"""
get domain

args:
    space::AbstractSpace
ret:
    AbstractDomain
"""
function domain end

"""
args:
    space::AbstractSpace{T,D}
ret:
    (x1, ..., xD,) # incl end points
"""
function points end

"""
args:
    space::AbstractSpace{T,D}
ret:
    Gauss quadratues in form
    ((z1,w1), ..., (zD,wD),)
"""
function quadratures end

"""
args:
    space::AbstractSpace{T,D}
ret:
    mass matrix as a SciMLOperator
"""
function mass_matrix end

"""
args:
    space::AbstractSpace{T,D}
ret:
    AbstractArray of size size(space)
"""
function local_numbering end

"""
args:
    space::AbstractSpace{T,D}
ret:
    AbstractArray of size size(space)
"""
function global_numbering end

"""
args:
    space::AbstractSpace{T,D}
    i::Integer
ret:
    ith basis function that can be evaluated
    anywhere in Space.domain
"""
function basis end

"""
get indices of boudnary nodes
"""
function boundary_nodes end

### interpolation
"""
Interpolate function values to to points.

args:
    points::vector of coordinates
    u::AbstractVector
    space::AbstractSpace{T,D}
ret:
    u interpolated to points
"""
function interp end

"""
Point-to-point interpolant between
spaces on same domain. used for
dealiasing

args:
    space1::AbstractSpace{T,D}
    space2::AbstractSpace{T,D}
ret:
    interpolation operator from
    space1 to space2
"""
function interpOp end

### modes

"""
modes
"""
function modes end

"""
size of modal basis
"""
function mode_size end

### transform

"""
Lazily transforms between physical and modal space.
"""
function transform end

"""
forward transform operator

args:
    - space::AbstractSpace
ret:
    - transform operator on space

to change transform operator to act on a different
AbstractVecOrMat type, call

`space = make_transform(space, input_vec)`
"""
function transformOp end

"""
zero high frequency modes

args:
    - space::AbstractSpace
    - frac::Number fraction of spectrun to keep
ret:
    - truncation operator on space
"""
function truncationOp end

"""
Form transform operator per new input vector for space

args:
    - space::AbstractSpace
    - u::AbstractVecOrMat

kwargs:
    - p: parameters
    - t: time

ret:
    - ftransform Forward transform wrapped
    `SciMLOperators.FunctionOperator`
"""
function form_transform end

"""
Set transform operator for space acting on vectors like u

args:
    - space::AbstractSpace
    - u::AbstractVecOrMat input prototype of size (N,) or (N,K)
      where N is the length(space). K will be the batch size.

ret:
    - space::AbstractSpace with transform operator that can act on `u`
"""
function make_transform(space::AbstractSpace,
                        u::Union{Nothing,AbstractArray} = nothing;
                        kwargs...)

    u = if u isa Union{AbstractVecOrMat,Nothing}
        u
    else # ND Array
        N = length(space)
        s = size(u)
        K = prod(s[2:end])

        @assert s[1] == N "Dimension mismatch"
        reshape(u, (N, K))
    end

    ftr = form_transform(space, u; kwargs...)
    @set! space.ftransform = ftr

    space
end
#
