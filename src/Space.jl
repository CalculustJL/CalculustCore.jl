#
# TODO
#   - rename "Space" to "Discretization" because
#       1. space is how you represent functions (with basis)
#       2. discretization is how you compute vector calculus operations

###
# AbstractSpace interface
###

"""
args:
    space::AbstractSpace{T,D}
ret:
    (x1, ..., xD,) # incl end points
"""
function get_grid end

"""
args:
    space::AbstractSpace{T,D}
ret:
    AbstractField{Integer, D}
"""
function get_global_numbering end

"""
get domain

args:
    space::AbstractSpace
ret:
    AbstractDomain
"""
function get_domain end

"""
args:
    space::AbstractSpace{T,D}
    i::Integer
ret:
    ith basis function that can be evaluated
    anywhere in Space.domain
"""
function get_basis end

"""
get number of points
"""
function get_numpoints end

### interpolation
"""
Interpolate function values to to points.

args:
    points::vector of coordinates
    u::AbstractField
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

Continuity isn't enforced across
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
    massOp: AbstractField -> AbstractField
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
    laplaceOp: AbstractField -> AbstractField
"""
function laplaceOp(space::AbstractSpace)
    D = gradOp(space)
    M = massOp(space)

    lapl = D' * [M] * D

    first(lapl)
end

"""
args:
    space::AbstractSpace{<:Number,D}
    vel...::AbstractField{<:Number,D}
    space_dealias
ret:
    advectionOp: AbstractField -> AbstractField

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
                     vel::AbstractField{<:Number,D}...
                    ) where{D}
    V = [DiagonalOp.(vel)...]

    Dx = gradOp(space)
    M  = massOp(space)

    advectOp = V' * [M] * Dx

    first(advectOp)
end

"""
Divergence
"""
function divergenceOp(space::AbstractSpace)
    D = gradOp(space)
    return reshape(G, 1, length(D))
end

### dealiased operators

function massOp(space1::AbstractSpace{<:Number,D},
                space2::AbstractSpace{<:Number,D};
                J = nothing,
               ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    M2 = massOp(space2)

    J12' ∘ M2 ∘ J12
end

function laplaceOp(space1::AbstractSpace{<:Number,D},
                   space2::AbstractSpace{<:Number,D};
                   J = nothing,
                  ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    M2 = massOp(space2)
    D1 = gradOp(space1)
    JD = [J12] .∘ D1

    laplOp = JD' * [M2] * JD

    first(laplOp)
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
                     vel::AbstractField{<:Number,D}...;
                     J = nothing,
                    ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    V1 = [DiagonalOp.(vel)...]
    V2 = J12 .* V1

    M2 = massOp(space2)
    D1 = gradOp(space1)

    advectOp = [J12]' * V2' * [M2] * [J12] * D1

    first(advectOp)
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
