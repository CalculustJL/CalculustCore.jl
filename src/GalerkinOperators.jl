#
###
# Gather-Scatter Operators - enforce continuity/ periodicity
###

Base.adjoint(A::AbstractGatherScatterOperator) = A

# TODO write GatherScatterOp that calls NNlib.gather, scatter
# implement mul!, *
struct GatherScatter{D} <: AbstractGatherScatterOperator{D}
    global_numbering
    implementation
end

"""
Q*Q'*u where Q: local -> global operator
"""
function DSS(u,l2g,g2l)

    Qu   = NNlib.scatter(+,u,l2g) # Q
    QQtu = NNlib.gather(Qu,g2l)   # Q'

    return v
end

"""
map from local vector to global vector
"""
function Qmatrix(n::Integer, periodic::Bool)
    Q = sparse(I,n, n-1)

    if periodic
        Q[end,1] = 1
    end

    Q
end

function GatherScatter(space::AbstractSpectralSpace{<:Number,D}) where{D}
    domain = get_domain(space)
    periodic = isperiodic(domain)
    npoints = get_numpoints(space)

    if !prod(periodic...)
        return IdentityOp{D}()
    end

    Qmats = Qmatrix.(npoints, periodic)

    Q = if D == 1
        MatrixOp(Qmats...)
    elseif D == 2
        TensorProductOp2D(Qmats...)
    elseif D == 3
        TensorProductOp3D(Qmats...)
    end

    QQt = Q * Q' # replace with call to NNlib gather-scatter
end

###
# Boundary Condition application
#
# Apply this boundary condition to that boundary tag
###

struct DirichletBC end

struct BoundaryCondition{T,D} <: AbstractBoundaryCondition{T,D}
    tags_to_type # dictionary
    type # dirichlet, neumann
    dirichlet_func! # (ub, space) -> mul!(ub, I, false)
#   neumann_func!
    mask # implementation
end

function BoundaryCondition(tags, space::AbstractSpace<:Number,2;
                           dirichlet_func! =nothing, neumann_func! = nothing)

    mask = generate_mask(tags, space)

    BoundaryCondition()
end

"""
 bc = (:Dirichlet,:Neumann,:Dirichlet,:Dirichlet) at (rmin, rmax, smin, smax)

 :Dirichlet = Dirichlet = zeros ∂Ω data\n
 :Neumann   = Neumann   = keeps ∂Ω data

 A periodic mesh overwrites 'D' to 'N' in direction of periodicity.

 To achieve inhomogeneous Dirichlet condition, apply the formulation
 u = ub + uh, where uh is homogeneous part, and ub is an arbitrary
 smooth function on Ω. Then, solve for uh
"""
function generate_mask(tags, space::AbstractSpace{<:Number,2})
    (nr, ns,) = space.npoints

    periodic = isperiodic(space.domain)

    Ix = sparse(I,nr,nr)
    Iy = sparse(I,ns,ns)

    ix = collect(1:(nr))
    iy = collect(1:(ns))

    if(bc[1] == :Dirichlet) ix = ix[2:end]   end
    if(bc[2] == :Dirichlet) ix = ix[1:end-1] end
    if(bc[3] == :Dirichlet) iy = iy[2:end]   end
    if(bc[4] == :Dirichlet) iy = iy[1:end-1] end

    if(periodic[1]) ix = collect(1:(nr)); end
    if(periodic[2]) iy = collect(1:(ns)); end

    Rx = Ix[ix,:]
    Ry = Iy[iy,:]

    M = diag(Rx'*Rx) * diag(Ry'*Ry)'
    M = Array(M) .== true

    return M
end

function applyBC!(u::AbstractField{<:Number,D}, bc::BoundaryCondition{<:Number,D}) where{D}

    return u
end

function applyBC!()
end
#
