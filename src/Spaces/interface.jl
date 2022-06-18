#
Base.eltype(::AbstractSpace{T}) where{T} = T

"""
Dimension of underlying domain
"""
PDEInterfaces.dims(::AbstractSpace{<:Any,D}) where{D} = D

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
plot of function over space
args:
    - u scalar field
    - space AbstractSpace
"""
Plots.plot

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
function grid end

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
