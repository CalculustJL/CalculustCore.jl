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
plot of function over space
args:
    - u scalar field
    - space AbstractSpace
"""
function Plots.plot(u::AbstractVector, space::AbstractSpace{<:Any,1}; kwargs...)

    (x,) = points(space)
    plt  = plot(x, u; kwargs...)
end

function Plots.plot(u::AbstractVector, space::AbstractSpace{<:Any,2}; a=30, b=30, kwargs...)

    npts = size(space)
    (x,y) = points(space)

    u = _reshape(u, npts)
    x = _reshape(x, npts)
    y = _reshape(y, npts)

#   plt = plot(x, y, u, legend=false, c=:grays, camera=(a,b))
#   plt = plot!(x', y', u', legend=false, c=:grays, camera=(a,b))

    plt = Plots.heatmap(u; kwargs...)

    plt
end

function Plots.animate(u::AbstractMatrix, space::AbstractSpace{<:Any,1}; kwargs...)
    ylims = begin
        u = sol.u[1]
        mi = minimum(u)
        ma = maximum(u)
        buf = (ma-mi)/5
        (mi-buf, ma+buf)
    end
    anim = @animate for i=1:size(u, 2)
        plt = plot(u[:,i], space; ylims=ylims, kwargs...)
    end
end

function Plots.animate(u::AbstractMatrix, space::AbstractSpace{<:Any,2}; kwargs...)
    clim = begin
        mi = minimum(u)
        ma = maximum(u)
        (mi, ma)
    end
    anim = @animate for i=1:size(u, 2)
        plt = plot(u[:,i], space; clim=clim, kwargs...)
    end
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
    - frac::Number truncation fraction
ret:
    - truncation operator on space
"""
function truncationOp end

"""
Form transform operator per new input vector for space

args:
    - u::AbstractVecOrMat
    - space::AbstractSpace
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
                        u::AbstractVecOrMat{T} = first(points(space));
                        p=nothing,
                        t=zero(T),
                       ) where{T}
    ftr = form_transform(u, space, p=p, t=t)
    @set! space.ftransform = ftr

    space
end
#
