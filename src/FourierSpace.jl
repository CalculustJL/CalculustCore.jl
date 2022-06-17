#
# look at FourierFlows.jl for ref
struct FourierSpace{
                    T,
                    D,
                    Tdom<:AbstractDomain{<:Number,D},
                    Tgrid,
                    Tmass,
                    Tderiv,
                    Tloc,
                   }
    """ Domain """
    domain::Tdom
    """ size tuple """
    npoints::Tpts
    """ grid points """
    grid::Tgrid
    """ mass matrix """
    mass_matrix::Tmass
    """ forward transform `mul!(û, T , u)` """
    Tr::Ttr
    """ inverse transform `mul!(u, T , û)` """
    iTr::Titr
end

function FourierSpace(npts...;
                      domain=reference_box(length(npts))
                     )

    FourierSpace(
                 domain, npoints, grid,
                )
end

struct Collocation <: AbstractDiscretization end
struct Galerkin    <: AbstractDiscretization end

function gradOp(space::FourierSpace)
    D = dims(space)
end

function laplaceOp(space::FourierSpace, ::Galerkin)
end

function laplaceOp(space::FourierSpace, ::Collocation)
end

### interpolation operators

function interpOp(space1::FourierSpace, space2::FourierSpace)
end

#
