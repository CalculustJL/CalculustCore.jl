#
import FFTW: plan_rfft, plan_irfft

# TODO allow for evolution in either transformed space, or physical space

# look at FourierFlows.jl for ref
struct FourierSpace{
                    T,
                    D,
                    Tdom<:AbstractDomain{<:Any,D},
                    Tgrid,
                    Tmass,
                    Ttr,
                    Titr,
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
    transform::Ttr
    """ inverse transform `mul!(u, T , û)` """
    itransform::Titr
end

function FourierSpace(n::Integer;
                      domain::AbstractDomain{<:Any,1}=reference_box(1),
                      T=Float64,
                     )

    if domain isa IntervalDomain
        domain = BoxDomain(domain)
    elseif !(domain isa BoxDomain)
        @error "Trigonometric polynomials work with logically rectangular domains"
    end

    domain = PeriodicBox(((-π, π),))
    #""" reset deformation to map from [-π,π]^D """
    #ref_domain = reference_box(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO

    domain = T(domain)
    npoints = (n,)
    grid = 
    mass_matrix = 
    transform = 
    itransform = 

    FourierSpace(
                 domain, npoints, grid,
                 mass_matrix, transform, itransform,
                )

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
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
