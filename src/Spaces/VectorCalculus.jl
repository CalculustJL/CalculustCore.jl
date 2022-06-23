
###
# Vector Calculus Operators
###

"""
Mass Operator - ∫

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
Gradient Operator - ∇

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
Hessian Operator - ∇²
"""
function hessianOp end
hessianOp(space::AbstractSpace, discr::AbstractDiscretization) = hessianOp(space)

"""
Laplace Operator - Δ

args:
    space::AbstractSpace{T,D}
    space_dealias
ret:
    laplaceOp: AbstractVector -> AbstractVector
"""
function laplaceOp end

"""
Diffusion operator - νΔ
"""
function diffusionOp(ν::Number, args...)
    ν * laplaceOp(args...) # TODO - ν update behaviour
end

"""
Diffusion operator - ∇⋅(ν∇⋅)
"""
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
    JDDt = _transp(JDD, discr)

    - JDDt * MMν2 * JDD
end

"""
Advection Operator - v⃗⋅∇

args:
    vel...::AbstractVector
    space::AbstractSpace{<:Any,D}
    space_dealias (optional)
ret:
    advectionOp: AbstractVector -> AbstractVector
"""
function advectionOp end

"""
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
                     discr::AbstractDiscretization;
                     vel_update_funcs=nothing,
                    ) where{D}

    VV = []
    for i=1:D
        vel_update_func = if vel_update_funcs isa Nothing
            DEFAULT_UPDATE_FUNC
        else
            vel_update_funcs[i]
        end

        function update_func!(A, u, p, t)
            vel_update_func(A.diag, u, p, t)
            A
        end
        V = MatrixOperator(Diagonal(vel[i]); update_func=update_func!)
        push!(VV, V)
    end

    DD = gradOp(space, discr)
    M  = massOp(space, discr)
    MM = Diagonal([M for i=1:D])

    VV' * MM * DD
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
Divergence Operator - ∇⋅
"""
function divergenceOp end

"""
Add forcing
"""
function forcing end

function forcing(f::AbstractVector, discr::AbstractDiscretization)
    M = massOp(space, discr)
    M * f
end
#
