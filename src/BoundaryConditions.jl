#
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
    x = get_grid(space) |> first

    antimasks = []
    for i=1:num_boundaries(domain)
        idx = indices[i]
        M = similar(x, Bool) .* false
        if prod(length.(idx)) == 1
            M[idx...] = true
        else
            M[idx...] .= true
        end
        push!(antimasks, M)
    end

    DiagonalOp.(antimasks)
end

function dirichlet_mask(space, domain, indices, bc_dict)
    tags = boundary_tags(domain)
    x    = get_grid(space) |> first
    M    = similar(x, Bool) .* false .+ true

    for i=1:num_boundaries(domain)
        tag = boundary_tag(domain, i)
        bc  = bc_dict[tag]

        if bc isa DirichletBC
            idx = indices[i]
            if prod(length.(idx)) == 1
                M[idx...] = false
            else
                M[idx...] .= false
            end
        end
    end

    DiagonalOp(M)
end

struct BoundaryCondition{T,D,Tbcs,Tamasks,Tmask,Tamask,
                         Tspace<:AbstractSpace{T,D},
                        } <: AbstractBoundaryCondition{T,D}
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
    amask_dir = IdentityOp(space) - mask_dir

    BoundaryCondition(bc_dict, antimasks, mask_dir, amask_dir, space)
end
#
