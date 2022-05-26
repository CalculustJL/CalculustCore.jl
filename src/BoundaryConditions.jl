#
###
# BC Types
###

"""
u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct DirichletBC{F}
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

"""
 A periodic mesh overwrites 'D' to 'N' in direction of periodicity.

 To achieve inhomogeneous Dirichlet condition, apply the formulation
 u = ub + uh, where uh is homogeneous part, and ub is an arbitrary
 smooth function on Ω. Then, solve for uh
"""

function dirichlet_mask(space::AbstractSpace{<:Number,D}, bcs) where{D}
    domain  = get_domain(space)
    indices = boundary_nodes(space)

    grid = get_grid(space)
    M = @. false * similar(grid[1], Bool) + true

    for i=1:2D
        tag = boundary_tag(domain, i)
        bc  = bcs[tag]
        idx = indices[i]

        # ignore BC if periodic
        if iperiodic(domain, i÷2+1)
            break
        end

        # mask Dirichlet BC
        if bc isa DirichletBC
            M[idx] = false
        end
    end

    DiagonalOp(M)
end

struct BoundaryCondition{T,D} <: AbstractBoundaryCondition{T,D}
    """Dict(Domain_bdry_tag => BCType) """
    bcs
    """Vector(Domain_Bdry_Tag => Node indices) """
    indices
    """Diagonal Mask operator hiding Dirichlet boundaries"""
    M
    """Function Space"""
    space
end

function BoundaryCondition(bcs, space::AbstractSpace{<:Number,2})
    BoundaryCondition()
end

function applyBC!(u::AbstractField{<:Number,D}, bc::BoundaryCondition{<:Number,D}) where{D}
    @unpack bcs, indices, M, space = bc

    u = M * u
    return u
end

function applyBC(u::AbstractField{<:Number,D}, bc::BoundaryCondition{<:Number,D}) where{D}
    return u
end
#
