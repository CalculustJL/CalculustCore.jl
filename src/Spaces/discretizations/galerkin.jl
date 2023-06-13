#
###
# Galerkin discretization
###

"""
$TYPEDEF

The `Galerkin` discretization is a weighted residual method where the trial
space is equal to the test space.
"""
struct Galerkin <: AbstractDiscretization end

"""
$SIGNATURES

```
[Dx  --> [Dx' Dy']
 Dy]
```
"""
_transp(v::AbstractVector, ::Galerkin) = adjoint(v)

# operators

GALERKIN_DEALIAS_MSG = """
    `Galerkin` discretizations suppport dealiased implementation of vector
    calculus operations where integration is performed in  a(usually higher
    order) space, `Vd`. This is referred to as over-integration. If `Vd` is not
    provided, integration is done in `V`. `J` is the grid-to-grid interpolation
    operator from `V` to `Vd`. If `J` is not provided, it is recomputed. If
    `Vd == V`, then `J` becomes `IdentityOperator(V)`.
    """

"""
$SIGNATURES

$GALERKIN_DEALIAS_MSG
"""
function massOp(V::AbstractSpace{<:Any, D},
                discr::Galerkin,
                Vd::AbstractSpace{<:Any, D},
                J = nothing) where {D}

    J = isnothing(J) ? interpOp(V, Vd) : J
    # Jt = _transp(J, discr)

    M2 = massOp(Vd, discr)

    J' * M2 * J
end

"""
$SIGNATURES

For a `Galerkin` discretization, the entries of the negative laplacian are
given by

``
[A]_{ij} = \\int_\\Omega \\nabla\\phi_i(x) \\cdot \\nabla\\phi_j(x) dx
``

for v,u in H¹₀(Ω)

(v,-∇² u) = (vx,ux) + (vy,uy)\n
         := a(v,u)\n
"""
function laplaceOp(V::AbstractSpace{T, D}, discr::Galerkin) where{T, D}

    M  = massOp(V, discr)
    DD = gradientOp(V, discr)
    MM = Diagonal([M  for _ in 1:D])

    DD' * MM * DD
end

"""
$SIGNATURES

$GALERKIN_DEALIAS_MSG
"""
function laplaceOp(V::AbstractSpace{<:Any, D},
                   discr::Galerkin,
                   Vd::AbstractSpace{<:Any, D},
                   J = nothing) where{D}

    J = isnothing(J) ? interpOp(V, Vd) : J

    Md = massOp(Vd, discr)
    DD = gradientOp(V, discr)

    JMJ = J' * M * J
    MM = Diagonal([JMJ  for _ in 1:D])

    DD' * MM * DD
end

"""
    biharmonicOp(V::AbstractSpace{<:Any, D}, discr::Galerkin, Vd::AbstractSpace{<:Any, D}, J = nothing) where{D}

$GALERKIN_DEALIAS_MSG
"""
biharmonicOp

"""
$SIGNATURES

$GALERKIN_DEALIAS_MSG
"""
function diffusionOp(ν::AbstractVecOrMat,
                     V::AbstractSpace{<:Any, D},
                     discr::Galerkin,
                     Vd::AbstractSpace{<:Any, D},
                     J = nothing;
                     ν_update_func = DEFAULT_UPDATE_FUNC) where{D}

    @error "TODO: dealiased galerkin diffusion not implemented"

    J = isnothing(J) ? interpOp(V, Vd) : J

    ν = DiagonalOperator(ν; update_func = ν_update_func)
    Jν = J * ν # <-- want to put J * ν in diagonal matrix. but need update
    # use function operator?

    M2 = massOp(Vd, discr)
    Mν2 = Jν * M2
    MMν2 = Diagonal([Mν2 for i in 1:D])

    DD = gradientOp(V, discr)
    JDD = J .* DD
    JDDt = _transp(JDD, discr)

    JDDt * MMν2 * JDD
end

"""
$SIGNATURES

$GALERKIN_DEALIAS_MSG

ux,uy, ∇T are interpolated to a grid with higher polynomial order
for dealiasing (over-integration) so we don't commit any "variational crimes"
"""
function advectionOp(vels::NTuple{D},
                     V::AbstractSpace{<:Any, D},
                     discr::Galerkin,
                     Vd::AbstractSpace{<:Any, D};
                     J = nothing;
                     vel_update_funcs = nothing) where {D}

    @error "TODO dealiased galerkin advection has not been implemented"

    J12 = J !== nothing ? J : interpOp(V1, V2)

    VV1 = [DiagonalOperator.(vels)...]
    VV2 = J12 .* VV1
    M2 = massOp(V2, discr)
    MM2 = Diagonal([M for i in 1:D])
    DD1 = gradientOp(V1)

    VV2' * MM2 * (J12 .* DD1)
end
#
