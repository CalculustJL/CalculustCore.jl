#
struct BoundaryCondition{T,
                         Tdict,
                         Tamasks,
                         Tmask,
                         Tamask,
                         Tspace<:AbstractSpace{T},
                         Tdiscr<:AbstractDiscretization,
                        } <: AbstractBoundaryCondition{T}
    """Dict(Domain_bdry_tag => BCType)"""
    bc_dict::Tdict
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
