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
    gradientOp(V::AbstractSpace{T, D}) where{T, D} -> [∂x1, ..., ∂xD]

Gradient Operator: ∇

Returns a size `D` `Vector` of operators. The `d`th operator corresponds to
the linear transformation from a function to its partial derivative in the
`d`th dimension at the grid points of `V`. Its `[ij]`th entry is equal to the 
value of the derivative of the `j`th basis function at the `i`th grid point.

``
[D]_{ij} = \\partial_{x}\\phi_{j}(x_i)
``
"""
gradientOp(V::AbstractSpace, ::AbstractDiscretization) = gradientOp(V)

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
hessianOp(V::AbstractSpace, ::AbstractDiscretization) = hessianOp(V)

function hessianOp(V::AbstractSpace)
    DD = gradientOp(V)
    DD_ = reshape(DD, (1, ndims(V)))

    DD * DD_
end

"""
    laplaceOp(V::AbstractSpace, discr::AbstractDiscretization)

Laplace Operator: -Δ

Returns the negative Laplace operator over `V` per discretization scheme
`discr`. As the Laplace operator (laplacian) is a negative-definitie operator,
we return the negative laplacian which is a positive-definite operator.
"""
function laplaceOp end

"""
    biharmonicOp(V::AbstractSpace, discr::AbstractDiscretization)

Biharmonic Operator: Δ²

Represents the Biharmonic opeator over `V` per discretization scheme `discr`.
"""
function biharmonicOp end

"""
    diffusionOp(ν::Number, V::AbstractSpace, discr::AbstractDiscretization)

Diffusion operator: -νΔ

Returns the negative Laplace Operator scaled by diffusion coefficient `ν`.
`ν` can be updated per `ν_update_func` which expects the same signature as
`SciMLOperators.ScalarOperator`. All positional arguments besides `ν` are
passed down to `laplaceOp` to form the negative lapalcian.
"""
function diffusionOp(ν::Number, args...; ν_update_func = DEFAULT_UPDATE_FUNC)

    ν = ScalarOperator(ν; update_func = ν_update_func)
    A = laplaceOp(args...)

    ν * A
end

"""
$SIGNATURES

Diffusion operator: -∇⋅(ν∇⋅)

Returns the diffusion operator where `ν` is the space-varying diffusion
coefficient in `V`. `ν` can be updated per `ν_update_func` which expets the
same signature as `SciMLOperators.DiagonalOperator`.

It is assumed that `ν` is a function in space `V`. However, if `V` is a
transformed space, then `ν` is a function in the physical space, i.e.
`transform(V)`.
"""
function diffusionOp(ν::AbstractVecOrMat,
                     V::AbstractSpace{<:Any, D},
                     discr::AbstractDiscretization;
                     ν_update_func = DEFAULT_UPDATE_FUNC) where{D}

    M = massOp(V, discr)
    ν = DiagonalOperator(ν; update_func = ν_update_func)

    Mν = M * ν
    MMν = Diagonal([Mν for _ in 1:D])

    DD = gradientOp(V, discr)

    # TODO - hack so ops play nice inside arrays
    DD = AbstractSciMLOperator[DD...]
    DDt = _transp(DD, discr)

    L = DDt * MMν * DD

    isa(L, AbstractArray) ? L[1] : L
end

"""
$SIGNATURES

Advection Operator: (v⃗⋅∇)⋅  

for v,u,T in H¹₀(Ω)

(v,(u⃗⋅∇)T) = (v,ux*∂xT + uy*∂yT)\n
           = v' *B*(ux*∂xT + uy*∂yT)\n

implemented as

(u⃗⋅∇)T = ux*∂xT + uy*∂yT

       = [ux uy] * [Dx] T
                   [Dx]
"""
function advectionOp(vels::NTuple{D},
                     V::AbstractSpace{<:Any, D},
                     discr::AbstractDiscretization;
                     vel_update_funcs = nothing) where {D}
    VV = _pair_update_funcs(vels, vel_update_funcs)
    VVt = _transp(VV, discr)

    DD = gradientOp(V, discr)
    M = massOp(V, discr)
    MM = Diagonal([M for i in 1:D])

    VVt * MM * DD
end

"""
    divergenceOp(V::AbstractSpace, discr::AbstractDiscretization)

Divergence Operator: ∇⋅
"""
function divergenceOp end

"""
$SIGNATURES

Added forcing as an operator.

F = forcingOp(f)
L = A + F

F(u) = 0 * u + M * f
L(u) = A * u + M * f
"""
function forcingOp(f::AbstractVecOrMat,
                   V::AbstractSpace,
                   discr::AbstractDiscretization;
                   f_update_func = DEFAULT_UPDATE_FUNC,
                   f_update_func! = DEFAULT_UPDATE_FUNC,
                  )
    Z = NullOperator(V)
    M = massOp(V, discr)

    AffineOperator(Z, M, f; update_func = f_update_func,
                   update_func! = f_update_func!)
end
#
