#

###
# abstract spaces interface
###

Base.eltype(::AbstractSpace{T}) where {T} = T

"""
$SIGNATURES

Dimension of underlying domain
"""
Base.ndims(::AbstractSpace{<:Any, D}) where {D} = D

"""
$SIGNATURES

get number of points
"""
Base.size(space::AbstractSpace{<:Any, D}, d) where {D} = size(space)[d] #TODO dims check

"""
$SIGNATURES

length of vector in space
"""
Base.length(space::AbstractSpace) = prod(size(space))

function SciMLOperators.IdentityOperator(space::AbstractSpace)
    N = length(space)
    SciMLOperators.IdentityOperator(N)
end

function SciMLOperators.NullOperator(space::AbstractSpace)
    N = length(space)
    SciMLOperators.NullOperator(N)
end

###
# interface to physical space
###

"""
    domain(V::AbstractSpace)

Get domain underlying `V`
"""
function domain end

"""
    points(V::AbstractSpace)


Get `NTuple{D}` of vectors corresonding to the grid
"""
function points end

"""
    quadratures(V::AbstractSpace)


Get underlying quadrature scheme
"""
function quadratures end

###
# boundary management interface
###

"""
    global_numbering(V::AbstractSpace)

Get `AbstractArray` of size `size(V)`
"""
function global_numbering end

"""
    boundary_nodes(V::AbstractSpace)

Get indices of boudnary nodes
"""
function boundary_nodes end

###
# interpolation interface
###
"""
    interp(points, u::AbstractVector, V::AbstractSpace)

Interpolate function described by `u` to points `points`.

"""
function interp end

"""
    interpOp(V1::AbstractSpace, V2::AbstractSpace)

Grid to grid interpolation operator from `V1` to `V2`. Both spaces must
share the same domain. Its `[ij]`th entry is the value of the `j`th basis
function of `V1` at the `i`th grid point of `V2`.

``
[J]_{ij} = \\phi_{j}(x_i)
``

If `V1 == V2`, we simply return the no-op `IdentityOperator(V1)`.
"""
function interpOp end

###
# spectral transform interface
###

"""
    modes(V::AbstractSpace)

Get `NTuple{D}` of vectors corresonding to the spectral grid
"""
function modes end

"""
    mode_size(V::AbstractSpace)

Get size of modal space. Equivalent to `size(transform(V))`
"""
function mode_size end

"""
    transform(V::AbstractSpace)

Lazily transforms between physical and modal space.
"""
function transform end

"""
    transformOp(V::AbstractSpace)

Get transform operator that takes vectors from `V` to `transform(V)`.
To change transform operator to act on a different `AbstractVecOrMat`
subtype, call

    V = make_transform(V, input_prototype)
"""
function transformOp end

"""
    truncationOp(V::AbstractSpace{T, D} , fracs::NTuple{D}) where{T, D}

Truncation operator removes high-frequency modes in an input vector.
`fracs[d]` is the proprtion of the spectrum to be preserved in `d`th
dimension. In spectral space, the operator corresponds to diagonal scaling,
where the first `fracs[d]` of the spectrum (low frequency modes) is
multiplied by unity, and the remainder (high-frequency modes) with zeros.
"""
function truncationOp end

"""
    form_transform(V::AbstractSpace, input::AbstractVecOrMat; [kwargs...])

Form transform operator given by `transformOp(V)` that may be applied to the
provided `input` prototype array per keyword arguments `kwargs`. `kwargs` may
include `p`, the parameter set.
"""
function form_transform end

"""
$SIGNATURES

Returns a copy of `V` with transform operator (given by `transformOp(V)`)
that may be applied to the provided `input` prototype array per provided
keyword arguments `kwargs`.
"""
function make_transform(V::AbstractSpace,
                        input::Union{Nothing, AbstractVecOrMat} = nothing;
                        kwargs...)

    input = isnothing(input) ? V.points[1] : input

    if length(V) != size(input, 1)
        msg = """First dimension of input prototype, $input,
            (of size $(size(u))) must equal the length of $V (of size
            $(size(V)))."""
        throw(DimensionMismatch(msg))
    end

    F = form_transform(V, input; kwargs...)
    @set! V.ftransform = F

    V
end
#
