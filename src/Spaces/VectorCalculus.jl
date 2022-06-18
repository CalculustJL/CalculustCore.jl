
###
# Vector Calculus Operators
###

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
    MM = Diagonal([M for i=1:D])

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

    VV = [DiagonalOperator.(vel)...]
    DD = gradOp(space)
    M  = massOp(space)
    MM = Diagonal([M for i=1:D])

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
    MM2 = Diagonal([M2 for i=1:D])

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

    VV1 = [DiagonalOperator.(vel)...]
    VV2 = J12 .* VV1
    M2  = massOp(space2)
    MM2 = Diagonal([M for i=1:D])
    DD1 = gradOp(space1)

    VV2' * MM2 * (J12 .* DD1)
end
#
