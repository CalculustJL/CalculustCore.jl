#

###
# Boundary Condition Types
###

"""
u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct DirichletBC{F}
    f::F = (points...) -> zero(first(points))
end

"""
(n⋅∇)u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct NeumannBC{F}
    f::F = (points...) -> zero(first(points))
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

function boundary_antimasks(space::AbstractSpace, dom::AbstractDomain, indices)
    loc_num = local_numbering(space)

    antimasks = []
    for i=1:num_boundaries(dom)
        M = similar(loc_num, Bool) * false |> _vec
        idx = indices[i]
        set_val!(M, true, idx)
        push!(antimasks, M)
    end

    DiagonalOperator.(antimasks)
end

function dirichlet_mask(space::AbstractSpace, dom::AbstractDomain, indices, bc_dict)
    tags = boundary_tags(dom)
    x    = grid(space) |> first
    M    = similar(x, Bool) * false .+ true |> _vec

    for i=1:num_boundaries(dom)
        tag = boundary_tag(dom, i)
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
    dom       = domain(space)
    indices   = boundary_nodes(space)
    antimasks = boundary_antimasks(space, dom, indices)
    mask_dir  = dirichlet_mask(space, dom, indices, bc_dict)
    amask_dir = IdentityOperator{length(space)}() - mask_dir

    BoundaryCondition(bc_dict, antimasks, mask_dir, amask_dir, space)
end
#
