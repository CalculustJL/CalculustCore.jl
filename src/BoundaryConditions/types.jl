#
"""
$TYPEDEF

$FIELDS
"""
struct BoundaryCondition{T,
                         Tdict,
                         Tamasks,
                         Tmask,
                         Tamask,
                         Tspace <: AbstractSpace{T},
                         Tdiscr <: AbstractDiscretization
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

"""
$SIGNATURES

"""
function BoundaryCondition(bc_dict::Dict, V::AbstractSpace,
                           discr::AbstractDiscretization)
    dom = domain(V)
    indices = boundary_nodes(V)
    antimasks = boundary_antimasks(V, dom, indices)
    mask_dir = dirichlet_mask(V, dom, indices, bc_dict)
    amask_dir = IdentityOperator(length(V)) - mask_dir

    BoundaryCondition(bc_dict, antimasks, mask_dir, amask_dir, V, discr)
end
#
