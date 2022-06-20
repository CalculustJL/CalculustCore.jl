#
import FFTW: plan_fft, plan_ifft
import FFTW: plan_rfft, plan_irfft

# switch to CUDA.plan_fft with GPU arrays

# might be possible to keep grid as range for multidim cases as well
# if not, switch to multidim arrays

"""
with stuff in physical space. operators leverage transform
"""
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
    transforms::Ttr
    """ inverse transform `mul!(u, T , û)` """
    itransforms::Titr
end

function FourierSpace(n::Integer;
                      domain::AbstractDomain{<:Any,1}=fourier_box(1),
                      T=Float64,
                     )

    if domain isa IntervalDomain
        domain = BoxDomain(domain)
    elseif !(domain isa BoxDomain)
        @error "Trigonometric polynomials work with logically rectangular domains"
    end

    domain = FourierDomain(1)
    L = length(domain)
    #""" reset deformation to map from [-π,π]^D """
    #ref_domain = reference_box(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO

    dx = L / n
    x = range(start=-L/2, stop=L/2-dx, length=n)

    if T <: Real
        k   = (0:n÷2-1) * (2π/L)
        tr  = plan_rfft()
       itr = plan_irfft()
    else
        k   = (-n÷2:n÷2-1) * (2π/L)
        tr  = plan_fft()
        itr = plan_ifft()
    end

    tr, itr = begin
        op = 
        tr = FunctionOperator(
                             )

        itr = FunctionOperator(
                              )
    end

    domain = T(domain)
    npoints = (n,)
    grid = (x,)
    modes = (k,)
    mass_matrix = 
    transforms = (tr,)
    itransforms = (itr,)

    FourierSpace(
                 domain, npoints, grid,
                 mass_matrix, transform, itransform,
                )

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

###
# interface
###

Base.size(space::LagrangePolynomialSpace) = space.npoints
domain(space::FourierSpace) = space.domain
points(space::FourierSpace) = space.grid
function quadratures(space::FourierSpace)
    N = length(space)
    L = length(domain(space))
    x = points(space) |> first
    w = ones(N) * (2π/L)

    ((x, w),)
end
modes(space::FourierSpace) = space.modes

###
# vector calculus
###

function massOp(space::FourierSpace{<:Any,1}, ::Collocation)
end

function gradOp(space::FourierSpace{<:Any,1}) # ∇
    k = modes(space)
    ik = im * DiagonalOperator(k)

    itr * ik * tr
end

function hessianOp(space::FourierSpace{<:Any,1}) # ∇²
    k = modes(space)
    ik2 = -DiagonalOperator(k .* k)

    itr * ik2 * tr
end

function laplaceOp(space::FourierSpace{<:Any,1}, ::Collocation)
    hessianOp(space) |> sum
end

###
# interpolation operators
###

# low-pass filter via restriction/extension matrix
function interpOp(space1::FourierSpace{<:Any,1}, space2::FourierSpace{<:Any,1})
    M = size(space2) # output
    N = size(space1) # input

    J = sparse(I, (M,N)) |> MatrixOperator
end
#
