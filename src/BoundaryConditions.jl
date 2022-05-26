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

"""
Periodic Boundary Condition
"""
struct PeriodicBC{Ttag}
    tag::Ttag
end

###
# BC implementation
###

"""
 A periodic mesh overwrites :Dirichlet to :Neumann in direction of periodicity.
 Periodicity is implemented via the gather-scatter operator.

 To achieve inhomogeneous Dirichlet condition, apply the formulation
 u = ub + uh, where uh is homogeneous part, and ub is an arbitrary
 smooth function on Ω. Then, solve for uh
"""

function dirichlet_mask(domain, bcs, indices)
    #periodic = isperiodic(domain)
    tags = boundary_tags(domain)

    x = get_grid(space)[1]
    M = Bool.(zero(x)) .+ true

    for i=1:num_boundaries(domain)
        tag = boundary_tag(domain, i)
        bc  = bcs[tag]
        idx = indices[i]

        if bc isa DirichletBC
            M[idx] = false
        end
    end

    DiagonalOp(M)
end

struct BoundaryCondition{T,D,Tbcs,Tidx,
                         Tmask::AbstractOperator{Bool,D},
                         Tspace::AbstractSpace{<:Number,D},
                        } <: AbstractBoundaryCondition{T,D}
    """Dict(Domain_bdry_tag => BCType)"""
    bcs::Tbcs
    """Vector(Domain_Bdry_Tag => Node indices)"""
    indices::Tidx
    """Diagonal Mask operator hiding Dirichlet boundaries"""
    mask::Tmask
    """Function space"""
    space::Tsp
end

function BoundaryCondition(bcs::Dict, space::AbstractSpace)
    domain = get_domain(space)
    indices = boundary_nodes(space)
    mask = dirichlet_mask(domain, bcs, indices)
    BoundaryCondition(bcs, indices, mask)
end
#
