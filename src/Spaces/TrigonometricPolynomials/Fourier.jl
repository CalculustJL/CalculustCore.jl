#
import FFTW: plan_rfft, plan_irfft

# TODO allow for evolution in either transformed space, or physical space
struct TransformedSpace
    space
end

function transform(space::FourierSpace)
    TransformedSpace(space)
end

function transform(space::TransformedSpace)
    space.space
end

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

    if T <: Real
        tr  = plan_rfft()
        itr = plan_irfft()
    else
        tr  = plan_fft()
        itr = plan_ifft()
    end

    domain = T(domain)
    npoints = (n,)
    grid = (x,)
    modes = (k,)
    mass_matrix = 
    transform = (tr,)
    itransform = (itr,)

    FourierSpace(
                 domain, npoints, grid,
                 mass_matrix, transform, itransform,
                )

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

###
# vector calculus
###

function gradOp(space::FourierSpace)
    D = dims(space)
end

###
# interpolation operators
###

function interpOp(space1::FourierSpace, space2::FourierSpace)
    # low-pass filter via restriction/matrix
end
#
