#

###
# BC implementation
###

function boundary_antimasks(space::AbstractSpace, dom::AbstractDomain, indices)
    loc_num = local_numbering(space)

    antimasks = []
    for i=1:num_boundaries(dom)
        M = similar(loc_num, Bool) * false |> vec
        idx = indices[i]
        set_val!(M, true, idx)
        push!(antimasks, M)
    end

    DiagonalOperator.(antimasks)
end

function dirichlet_mask(space::AbstractSpace, dom::AbstractDomain, indices, bc_dict)
    tags = boundary_tags(dom)
    x    = points(space) |> first
    M    = similar(x, Bool) * false .+ true |> vec

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
                         Tdiscr<:AbstractDiscretization,
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
    """ Discretization """
    discr::Tdiscr
end

function BoundaryCondition(bc_dict::Dict,
                           space::AbstractSpace{<:Number,D},
                           discr::AbstractDiscretization) where{D}

    dom       = domain(space)
    indices   = boundary_nodes(space)
    antimasks = boundary_antimasks(space, dom, indices)
    mask_dir  = dirichlet_mask(space, dom, indices, bc_dict)
    amask_dir = IdentityOperator{length(space)}() - mask_dir

    BoundaryCondition(bc_dict, antimasks, mask_dir, amask_dir, space, discr)
end
#
