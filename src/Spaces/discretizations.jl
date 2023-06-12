#
"""
$TYPEDEF

The `Galerkin` discretization is a weighted residual method where the trial
space is equal to the test space.
"""
struct Galerkin <: AbstractDiscretization end

_transp(a, ::Galerkin) = adjoint(a)

"""
$SIGNATURES

For a `Galerkin` discretization, the Laplacian

``
[A]_{ij} = \\int_\\Omega \\nabla\\phi_i(x) \\cdot \\nabla\\phi_j(x) dx
``
"""
function laplaceOp(V::AbstractSpace, discr::Galerkin)
    D = ndims(V)

    M = massOp(V, discr)
    MM = Diagonal([M for i in 1:D])

    DD = gradientOp(V, discr)

    -DD' * MM * DD
end

"""
$SIGNATURES

Form the laplace opeator over `V1` with integration performed in (usually
higher-order space) `V2`.

For a `Galerkin` discretization,

for v,u in H¹₀(Ω)

(v,-∇² u) = (vx,ux) + (vy,uy)\n
         := a(v,u)\n
"""
function laplaceOp(V1::AbstractSpace{<:Any, D},
                   V2::AbstractSpace{<:Any, D},
                   discr::Galerkin;
                   J = nothing) where {D}

    J12 = J !== nothing ? J : interpOp(V1, V2)
    #J21 = _transp(J12) # or interpOp(V2, V1) # TODO

    M2 = massOp(V2, discr)
    MM2 = Diagonal([M2 for i in 1:D])

    DD1 = gradientOp(V1, discr)
    JDD = J12 .* DD1

    JDD' * MM2 * JDD
end

function diffusionOp(ν::AbstractVector, V::AbstractSpace, discr::Galerkin)
    D = ndims(V)
    ν = DiagonalOperator(ν)
    DD = gradientOp(V)
    M = massOp(V, discr)
    Mν = ν * M
    MMν = Diagonal([Mν for i in 1:D])

    DD' * MMν * DD
end

"""
$TYPEDEF

A `Collocation` discretization is when the strong form of the equation is
discretized as it is.
"""
struct Collocation <: AbstractDiscretization end

_transp(a, ::Collocation) = reshape(a, (1, length(a)))

"""
For a `Collocation` discretization, the mass operator returns
`IdentityOperator(V)`.

"""
massOp(V::AbstractSpace, ::Collocation) = IdentityOperator(V)

function laplaceOp(V::AbstractSpace, discr::Collocation)
    DD2 = hessianOp(V, discr)
    -tr(DD2) # trace
end

function biharmonicOp(V::AbstractSpace, discr::Collocation)
    A = laplaceOp(V, discr)

    A * A
end

function divergenceOp(V::AbstractSpace, ::Collocation)
    DD = gradientOp(V)

    sum(DD)
end
#
