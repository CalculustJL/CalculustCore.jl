#
###
# Vector Calculus Operators
###

"""
    massOp(V::AbstractSpace, discr::AbstractDiscretization)

Mass Operator: ∫⋅dx

Represents the mass matrix for discretization scheme `discr`.

For a `Galerkin` discretization, the `[ij]`th entry corresponds to the
inner product of the `i`th and `j`th basis functions.

``
[M]_{ij} = \\int_\\Omega \\phi_i(x)\\phi_j(x) dx = \\langle \\phi_i, \\phi_j \\rangle
``
"""
function massOp end

"""
$SIGNATURES


"""
function massOp(V1::AbstractSpace{<:Any, D},
                V2::AbstractSpace{<:Any, D},
                discr::AbstractDiscretization;
                J = nothing) where {D}
    @error "this method has not been implemented yet"
    J12 = J !== nothing ? J : interpOp(V1, V2)
    #J21 = _transp(J12) # or interpOp(V2, V1) # TODO

    M2 = massOp(V2, discr)

    J12 * M2 * J12
end

"""
    gradientOp(V::AbstractSpace{T, D}) where{T, D} -> [∂x1, ..., ∂xD]

Gradient Operator: ∇

Returns a size `D` `Vector` of operators. The `d`th operator corresponds to
the linear transformation from a function to its partial derivative in the
`d`th dimension at the grid points of `V`. Its `[ij]`th entry is equal to the 
value of the derivative of the `j`th basis function at the `i`th grid point.

``
[D]_{ij} = \\\partial_{x}phi_{j}(x_i)
``
"""
gradientOp(space::AbstractSpace, discr::AbstractDiscretization) = gradientOp(space)

"""
    hessianOp(V::AbstractSpace{T, D}) where{T, D} ->

``
[∂x1∂x1 ... ∂x1∂xD,
        ...
 ∂xD∂x1 ... ∂xD∂xD ]
``

Hessian Operator: ∇²

Returns a `D × D` `Matrix` of operators whose `[ij]`th entry represents the
transformation from a function to its second partial derivative in the
`[ij]`th directions.
"""
function hessianOp end
hessianOp(space::AbstractSpace, ::AbstractDiscretization) = hessianOp(space)

function hessianOp(V::AbstractSpace)
    DD = gradientOp(V)

    DD * transpose(DD)
end

"""
Laplace Operator: -Δ

args:
    space::AbstractSpace{T,D}
    space_dealias
ret:
    laplaceOp: AbstractVector -> AbstractVector
"""
function laplaceOp end

"""
Biharmonic Operator: Δ²

args:
    space::AbstractSpace{T,D}
    space_dealias (optional)
ret:
    biharmonicOp: AbstractVector -> AbstractVector
"""
function biharmonicOp end

"""
Diffusion operator: -νΔ
"""
function diffusionOp(ν::Number, args...; ν_update_func = DEFAULT_UPDATE_FUNC)
    ν = ScalarOperator(ν; update_func = ν_update_func)
    A = laplaceOp(args...)

    ν * A
end

"""
Diffusion operator: -∇⋅(ν∇⋅)
"""
function diffusionOp(ν::AbstractVecOrMat,
                     space1::AbstractSpace{<:Any, D},
                     space2::AbstractSpace{<:Any, D},
                     discr::AbstractDiscretization;
                     J = nothing) where {D}
    J12 = J !== nothing ? J : interpOp(space1, space2)
    #J21 = _transp(J12) # or interpOp(space2, space1) # TODO

    Jν = J * DiagonalOperator(ν)

    M2 = massOp(space2, discr)
    Mν2 = Jν * M2
    MMν2 = Diagonal([Mν2 for i in 1:D])

    DD = gradientOp(space1, discr)
    JDD = J .* DD
    JDDt = _transp(JDD, discr)

    JDDt * MMν2 * JDD
end

"""
Advection Operator: v⃗⋅∇

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

(u⃗⋅∇)T = ux*∂xT + uy*∂yT

       = [ux uy] * [Dx] T
                   [Dx]
"""
function advectionOp(vels::NTuple{D},
                     space::AbstractSpace{<:Any, D},
                     discr::AbstractDiscretization;
                     vel_update_funcs = nothing) where {D}
    VV = _pair_update_funcs(vels, vel_update_funcs)

    DD = gradientOp(space, discr)
    M = massOp(space, discr)
    MM = Diagonal([M for i in 1:D])

    VV' * MM * DD # TODO - transpose instead of adjoint VVt=_transp(VV, discr)
end

"""
ux,uy, ∇T are interpolated to
a grid with higher polynomial order
for dealiasing (over-integration)
so we don't commit any
"variational crimes"
"""
function advectionOp(vel::NTuple{D},
                     space1::AbstractSpace{<:Any, D},
                     space2::AbstractSpace{<:Any, D};
                     J = nothing) where {D}
    @error "this method has not been implemented yet"
    J12 = J !== nothing ? J : interpOp(space1, space2)
    #J21 = _transp(J12) # or interpOp(space2, space1) # TODO

    VV1 = [DiagonalOperator.(vel)...]
    VV2 = J12 .* VV1
    M2 = massOp(space2)
    MM2 = Diagonal([M for i in 1:D])
    DD1 = gradientOp(space1)

    VV2' * MM2 * (J12 .* DD1)
end

"""
Divergence Operator: ∇⋅
"""
function divergenceOp end

"""
Added forcing as an operator

F = forcingOp(f)

F(u) = u + M*f
"""
function forcingOp(f::AbstractVecOrMat,
                   space::AbstractSpace,
                   discr::AbstractDiscretization;
                   f_update_func = DEFAULT_UPDATE_FUNC,
                   f_update_func! = DEFAULT_UPDATE_FUNC,
                  )
    Z = NullOperator(space)
    M = massOp(space, discr)

    AffineOperator(Z, M, f; update_func = f_update_func,
                   update_func! = f_update_func!)
end
#
