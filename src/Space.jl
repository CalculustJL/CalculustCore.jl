#
# TODO
#   - rename "Space" to "Discretization" because
#       1. space is how you represent functions (with basis)
#       2. discretization is how you compute vector calculus operations

###
# AbstractSpace interface
###

import Base: length, summary, size

"""
get domain

args:
    space::AbstractSpace
ret:
    AbstractDomain
"""
function get_domain end # TODO rename to domain

"""
args:
    space::AbstractSpace{T,D}
ret:
    (x1, ..., xD,) # incl end points
"""
function get_grid end # TODO rename to grid

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
get number of points
"""
Base.size

"""
length of vector in space
"""
Base.length(space::AbstractSpace) = prod(size(space))

"""
plot of function over space
args:
    - u scalar field
    - space AbstractSpace
"""
Plots.plot

"""
get indices of boudnary nodes
"""
function boundary_nodes end

function Base.summary(io::IO, space::AbstractSpace{T,D}) where{T,D}
    type_color, no_color = SciMLBase.get_colorizers(io)
    print(io,
          type_color, nameof(typeof(space)),
          no_color," over domain ",
          type_color,typeof(get_domain(space)),
          no_color," with uType ",
          type_color,typeof(first(get_grid(space))),
          no_color
         )
end

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

### vector calculus ops
"""
Gradient Operator
Compute gradient of u∈H¹(Ω).

Continuity isn't necessarily enforced across
element boundaries for gradients

args:
    space::AbstractSpace{T,D}
ret:
    gradOp: u -> [dudx1, ..., dudxD]
"""
function gradOp end

"""
Mass Operator

Inner Product: <u, v := u' * M * v

args:
    space::AbstractSpace{T,D}
    space_dealias
ret:
    massOp: AbstractVector -> AbstractVector
"""
function massOp end

"""
Laplace Operator

for v,u in H¹₀(Ω)
get_
(v,-∇² u) = (vx,ux) + (vy,uy)\n
         := a(v,u)\n

args:
    space::AbstractSpace{T,D}
    space_dealias
ret:
    laplaceOp: AbstractVector -> AbstractVector
"""
function laplaceOp(space::AbstractSpace{<:Number,D}) where{D}
    DD = gradOp(space)
    M  = massOp(space)
    MM = Diagonal(AbstractSciMLOperator[M for i=1:D])

    -(DD' * MM * DD)
end

"""
args:
    space::AbstractSpace{<:Number,D}
    vel...::AbstractVector
    space_dealias
ret:
    advectionOp: AbstractVector -> AbstractVector

for v,u,T in H¹₀(Ω)

(v,(u⃗⋅∇)T) = (v,ux*∂xT + uy*∂yT)\n
           = v' *B*(ux*∂xT + uy*∂yT)\n

implemented as

R'R * QQ' * B * (ux*∂xT + uy*∂yT)

(u⃗⋅∇)T = ux*∂xT + uy*∂yT

       = [ux uy] * [Dx] T
                   [Dx]
"""
function advectionOp(space::AbstractSpace{<:Number,D},
                     vel::AbstractVector...) where{D}

    VV = AbstractSciMLOperator[DiagonalOperator.(vel)...]
    DD = gradOp(space)
    M  = massOp(space)
    MM = Diagonal(AbstractSciMLOperator[M for i=1:D])

    VV' * MM * DD
end

"""
Divergence
"""
function divergenceOp(space::AbstractSpace{<:Number,D}) where{D}
    Dx = gradOp(space)
    return _reshape(Dx, 1, D)
end

### dealiased operators

function massOp(space1::AbstractSpace{<:Number,D},
                space2::AbstractSpace{<:Number,D};
                J = nothing,
               ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    M2 = massOp(space2)

    J12' * M2 * J12
end

function laplaceOp(space1::AbstractSpace{<:Number,D},
                   space2::AbstractSpace{<:Number,D};
                   J = nothing,
                  ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    M2  = massOp(space2)
    MM2 = Diagonal(AbstractSciMLOperator[M2 for i=1:D])

    DD1 = gradOp(space1)
    JDD = J12 .* DD1

    -(JDD' * MM2 * JDD)
end

"""
ux,uy, ∇T are interpolated to
a grid with higher polynomial order
for dealiasing (over-integration)
so we don't commit any
"variational crimes"
"""
function advectionOp(space1::AbstractSpace{<:Number,D},
                     space2::AbstractSpace{<:Number,D},
                     vel::AbstractVector...;
                     J = nothing,
                    ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    VV1 = AbstractSciMLOperator[DiagonalOperator.(vel)...]
    VV2 = J12 .* VV1
    M2  = massOp(space2)
    MM2 = Diagonal(AbstractSciMLOperator[M for i=1:D])
    DD1 = gradOp(space1)

    VV2' * MM2 * (J12 .* DD1)
end

###
# AbstractSpectralSpace
###

###
# Tensor Product Spaces
###

#struct TensorProductSpace{T,D1+D2,
#                          Tspace1<:AbstractSpace{<:Number, D1},
#                          Tspace2<:AbstractSpace{<:Number, D2},
#                         } <: AbstractTensorProductSpace{T,D1+D2}
#    space1::Tspace1
#    space2::Tspace2
#end
#
