#
import FFTW: plan_fft, plan_ifft, fftfreq
import FFTW: plan_rfft, plan_irfft, rfftfreq

# switch to CUDA.plan_fft with GPU arrays

"""
with stuff in physical space. operators leverage transform
"""
struct FourierSpace{
                    T,
                    D,
                    Tdom<:AbstractDomain{T,D},
                    Tpts,
                    Tgrid,
                    Tmodes,
                    Tmass,
                    Ttr,
                   }
    """ Domain """
    domain::Tdom
    """ size tuple """
    npoints::Tpts
    """ grid points """
    grid::Tgrid
    """ modes """
    modes::Tmodes
    """ mass matrix """
    mass_matrix::Tmass
    """ forward transform `mul!(û, T , u)` """
    transforms::Ttr
end

function FourierSpace(n::Integer;
                      domain::AbstractDomain{<:Any,1}=FourierDomain(1),
                      T=Float64,
                     )

    if domain isa IntervalDomain
        domain = BoxDomain(domain)
    elseif !(domain isa BoxDomain)
        @error "Trigonometric polynomials work with logically rectangular domains"
    end

    domain = FourierDomain(1)
    L = 2π #size(domain)
    #""" reset deformation to map from [-π,π]^D """
    #ref_domain = reference_box(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO

    ComplexT = if T isa Type{Float16}
        ComplexF16
    elseif T isa Type{Float32}
        ComplexF32
    else
        ComplexF64
    end

    dx = L / n
    x  = range(start=-L/2, stop=L/2-dx, length=n) |> Array

    if T <: Real
        k   = rfftfreq(n, 2π*n/L) |> Array
        tr  = plan_rfft(x)
        itr = plan_irfft(im*k, n)
    else
        k   = fftfreq(n, 2π*n/L) |> Array
        tr  = plan_fft(x)
        itr = plan_ifft(k, n)
    end

    transform = FunctionOperator(
                                 (du,u,p,t) -> mul!(du, tr, u);
                                 isinplace=true,
                                 T=ComplexT,
                                 size=(length(k),n),

#                                op_adjoint = 
                                 op_inverse = (du,u,p,t) -> ldiv!(du, tr, u)
                                )

    domain = T(domain)
    npoints = (n,)
    grid = (x,)
    modes = k #(k,)
    mass_matrix = ones(T, n) * (2π/L)
    transforms = transform#(transform,)

    space = FourierSpace(
                         domain, npoints, grid, modes,
                         mass_matrix, transforms,
                        )

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

###
# interface
###

Base.size(space::FourierSpace) = space.npoints
domain(space::FourierSpace) = space.domain
points(space::FourierSpace) = space.grid
function quadratures(space::FourierSpace{<:Any,1})
    x = points(space) |> first
    w = mass_matrix(space)

    ((x, w),)
end
mass_matrix(space::FourierSpace) = space.mass_matrix
modes(space::FourierSpace) = space.modes
transforms(space::FourierSpace) = space.transforms

## TODO - local system <-> global system
## global system for computation
## local system for plotting
# local_numbering(space::FourierSpace)
# global_numbering(space::FourierSpace)

###
# vector calculus
###

function massOp(space::FourierSpace{<:Any,1}, ::Collocation)
end

function gradOp(space::FourierSpace{<:Any,1}) # ∇
    tr = transforms(space)

    k  = modes(space)
    ik = im * DiagonalOperator(k)

    [
     tr \ ik * tr,
    ]
end

function hessianOp(space::FourierSpace{<:Any,1}) # ∇²
    tr = transforms(space)

    k   = modes(space)
    ik2 = -DiagonalOperator(@. k * k)

    [
     tr \ ik2 * tr,
    ]
end

function laplaceOp(space::FourierSpace{<:Any,1}, ::Collocation)
    hessianOp(space) |> sum
end

###
# interpolation operators
###

function interpOp(space1::FourierSpace{<:Any,1}, space2::FourierSpace{<:Any,1})
    tr1 = transforms(space1)

    M = length.(modes(space2)[1]) # output
    N = length.(modes(space1)[1]) # input

    J = sparse(I, (M,N)) |> MatrixOperator

    tr1 \ J * tr1
end
#
