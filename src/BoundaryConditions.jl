#
###
# BC Types
###

"""
u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct DirichletBC{F}
    f::F = grid = zero(grid[1])
end

"""
(n⋅∇)u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct NeumannBC{F}
    f::F = grid = zero(grid[1])
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

function boundary_antimasks(domain, indices)
    x = get_grid(space) |> first

    antimasks = []
    for i=1:num_boundaries(domain)
        idx = indices[i]
        M = similar(x, Bool) .* false
        M[idx] = true

        push!(antimasks, M)
    end

    DiagonalOp.(antimasks)
end

function dirichlet_mask(domain, bcs, indices)
    tags = boundary_tags(domain)
    x    = get_grid(space) |> first
    mask = Bool.(zero(x)) .+ true

    for i=1:num_boundaries(domain)
        tag = boundary_tag(domain, i)
        bc  = bcs[tag]

        if bc isa DirichletBC
            idx = indices[i]
            mask[idx] = false
        end
    end

    DiagonalOp(mask)
end

struct BoundaryCondition{T,D,Tbcs,Tamasks,
                         Tmask<:AbstractOperator{Bool,D},
                         Tamask<:AbstractOperator{Bool,D},
                         Tspace<:AbstractSpace{T,D},
                        } <: AbstractBoundaryCondition{T,D}
    """Dict(Domain_bdry_tag => BCType)"""
    bcs::Tbcs
    """Vector(boundary_antimasks); antimask = id - mask"""
    antimasks::Tamasks
    """Diagonal Mask operator hiding Dirichlet boundaries"""
    mask_dir::Tmask
    """antimask for dirichlet BC"""
    amask_dir::Tamask
    """Function space"""
    space::Tspace
end

function BoundaryCondition(bcs::Dict, space::AbstractSpace)
    domain    = get_domain(space)
    indices   = boundary_nodes(space)
    antimasks = boundary_antimasks(domain, indices)
    mask_dir  = dirichlet_mask(domain, bcs, indices)
    amask_dir = IdentityOp{D}() - mask_dir

    BoundaryCondition(bcs, bc_antimasks, mask_dir, amask_dir, space)
end
#
