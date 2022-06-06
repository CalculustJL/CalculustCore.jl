#
abstract type AbstractBoundaryCondition{T} end

###
# BC Types
###

"""
u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct DirichletBC{F}
    f::F = (grid...) -> zero(first(grid))
end

"""
(n⋅∇)u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct NeumannBC{F}
    f::F = (grid...) -> zero(first(grid))
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

function boundary_antimasks(space, domain, indices)
    loc_num = local_numbering(space)

    antimasks = []
    for i=1:num_boundaries(domain)
        M = similar(loc_num, Bool) * false |> _vec
        idx = indices[i]
        set_val!(M, true, idx)
        push!(antimasks, M)
    end

    DiagonalOperator.(antimasks)
end

function dirichlet_mask(space, domain, indices, bc_dict)
    tags = boundary_tags(domain)
    x    = get_grid(space) |> first
    M    = similar(x, Bool) * false .+ true |> _vec

    for i=1:num_boundaries(domain)
        tag = boundary_tag(domain, i)
        bc  = bc_dict[tag]

        if bc isa DirichletBC
            idx = indices[i]
            set_val!(M, false, idx)
        end
    end

    DiagonalOperator(M)
end

struct BoundaryCondition{T,Tbcs,Tamasks,Tmask,Tamask,
                         Tspace<:AbstractSpace{T},
                        } <: AbstractBoundaryCondition{T}
    """Dict(Domain_bdry_tag => BCType)"""
    bc_dict::Tbcs
    """Vector(boundary_antimasks); antimask = id - mask"""
    antimasks::Tamasks
    """Diagonal Mask operator hiding Dirichlet boundaries"""
    mask_dir::Tmask
    """antimask for dirichlet BC"""
    amask_dir::Tamask
    """Function space"""
    space::Tspace
end

function BoundaryCondition(bc_dict::Dict, space::AbstractSpace{<:Number,D}) where{D}
    domain    = get_domain(space)
    indices   = boundary_nodes(space)
    antimasks = boundary_antimasks(space, domain, indices)
    mask_dir  = dirichlet_mask(space, domain, indices, bc_dict)
    amask_dir = IdentityOperator{length(space)}() - mask_dir

    BoundaryCondition(bc_dict, antimasks, mask_dir, amask_dir, space)
end
#
