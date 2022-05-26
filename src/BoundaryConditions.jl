#
###
# BC Types
###

"""
u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwedf struct DirichletBC{F}
    f::F = zero
end

"""
(n⋅∇)u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct NeumannBC{F}
    f::F = zero
end

"""
    f1(x)u(x) + f2(x)(n⋅∇)u(x) = f3(x), x ∈ ∂Ω
"""
struct RobinBC{F1,F2,F3}
    f1::F1
    f2::F2
    f3::F3
end

###
# Boundary tag to index
###

"""
"""
struct Boundary_to_index_map
end

###
# BC Application
###

"""
 Apply this boundary condition to that boundary tag

mark boundary nodes with proper tag

basically

BCType --> Tag --> Node indices

Mapping 1: "bcs"
    Domain_BC_tag => BC_type (dirichlet, neumann, etc)

Mapping 2
    Domain_BC_tag => Node indices
"""

struct BCMapping{T,D} <: AbstractBoundaryCondition{T,D}
    """ Dict(Domain_Bdry_Tag => BoundaryCondition) """
    bcs
    """ Dict(Domain_Bdry_Tag => Node indices) """
    idx
    domain
    space
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
