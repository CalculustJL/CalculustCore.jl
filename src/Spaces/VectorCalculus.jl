
###
# Vector Calculus Operators
###

"""
Gradient Operator
Compute gradient of u∈H¹(Ω).

Continuity isn't necessarily enforced across
element boundaries for gradients

args:
    space::AbstractSpace
    discr::AbstractDiscretization (optional)
ret:
    gradOp: u -> [dudx1, ..., dudxD]
"""
function gradOp end
gradOp(space::AbstractSpace, discr::AbstractDiscretization) = gradOp(space)

"""
Mass Operator

Inner Product: <u, v := u' * M * v

args:
    space::AbstractSpace
    discr::AbstractDiscretization (optional)
    space_dealias
ret:
    massOp: AbstractVector -> AbstractVector
"""
function massOp end

function massOp(space1::AbstractSpace{<:Any,D},
                space2::AbstractSpace{<:Any,D},
                discr::AbstractDiscretization;
                J = nothing,
               ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)
    #J21 = _transp(J12) # or interpOp(space2, space1) # TODO

    M2 = massOp(space2)

    J12 * M2 * J12
end

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
function laplaceOp(space::AbstractSpace, discr::AbstractDiscretization)
    D = dims(space)

    M  = massOp(space, space)
    MM = Diagonal([M for i=1:D])

    DD = gradOp(space, discr)
    DDt = _transp(DD, disscr)

    - DDt * MM * DD
end

function laplaceOp(space1::AbstractSpace{<:Any,D},
                   space2::AbstractSpace{<:Any,D},
                   discr::AbstractDiscretization;
                   J = nothing,
                  ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)
    #J21 = _transp(J12) # or interpOp(space2, space1) # TODO

    M2  = massOp(space2, discr)
    MM2 = Diagonal([M2 for i=1:D])

    DD1  = gradOp(space1, discr)
    JDD  = J12 .* DD1
    JDDt = _transp(JDD, discr)

    - JDDt * MM2 * JDD
end

"""
Diffusion operator
"""
function diffusionOp(ν::Number, args...)
    ν * laplaceOp(args...)
end

function diffusionOp(ν::AbstractVector, space::AbstractSpace, discr::AbstractDiscretization)
    D = dims(space)
    ν = DiagonalOperator(ν)
    DD = gradOp(space)
    M  = massOp(space)
    Mν = ν * M
    MMν = Diagonal([Mν for i=1:D])
    DDt = _transp(DD, discr)

    - DDt * MMν * DD
end

function diffusionOp(ν::AbstractVector,
                     space1::AbstractSpace{<:Any,D},
                     space2::AbstractSpace{<:Any,D},
                     discr::AbstractDiscretization;
                     J = nothing,
                    ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)
    #J21 = _transp(J12) # or interpOp(space2, space1) # TODO

    Jν = J * DiagonalOperator(ν)

    M2  = massOp(space2, discr)
    Mν2 = Jv * M2
    MMν2 = Diagonal([Mν2 for i=1:D])

    DD  = gradOp(space, discr)
    JDD = J .* DD

    - JDDt * MMν2 * JDD
end

"""
args:
    space::AbstractSpace{<:Any,D}
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
function advectionOp(vel::NTuple{D},
                     space::AbstractSpace{<:Any,D},
                     discr::AbstractDiscretization
                    ) where{D}

    VV = [DiagonalOperator.(vel)...]
    DD = gradOp(space, discr)
    M  = massOp(space, discr)
    MM = Diagonal([M for i=1:D])

    VVt = _transp(VV, discr)

    VVt * MM * DD
end

"""
ux,uy, ∇T are interpolated to
a grid with higher polynomial order
for dealiasing (over-integration)
so we don't commit any
"variational crimes"
"""
function advectionOp(space1::AbstractSpace{<:Any,D},
                     space2::AbstractSpace{<:Any,D},
                     vel::AbstractVector...;
                     J = nothing,
                    ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)
    #J21 = _transp(J12) # or interpOp(space2, space1) # TODO

    VV1 = [DiagonalOperator.(vel)...]
    VV2 = J12 .* VV1
    M2  = massOp(space2)
    MM2 = Diagonal([M for i=1:D])
    DD1 = gradOp(space1)

    VV2' * MM2 * (J12 .* DD1)
end

"""
Divergence Operator
"""
function divergenceOp(space::AbstractSpace)
    D = dims(space)
    Dx = gradOp(space)
    return _reshape(Dx, (1, D))
end
#
