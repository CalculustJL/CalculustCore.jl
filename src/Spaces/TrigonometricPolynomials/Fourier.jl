#
import FFTW: plan_fft, plan_ifft, fftfreq
import FFTW: plan_rfft, plan_irfft, rfftfreq

"""
with stuff in physical space. operators use FFT
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
                   } <: AbstractSpace{T,D}
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
    ftransform::Ttr
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
    k  = rfftfreq(n, 2π*n/L) |> Array

    domain = T(domain)
    npoints = (n,)
    grid = (x,)
    modes = (k,)
    mass_matrix = ones(T, n) * (2π/L)
    ftransform = nothing

    space = FourierSpace(
                         domain, npoints, grid, modes,
                         mass_matrix, ftransform,
                        )

    space = make_transform(space, x)

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

function form_transform(u::AbstractVecOrMat, space::FourierSpace{T,D}) where{T,D}

    ssp = size(space)
    N   = length(space)

    @assert size(u, 1) == N "size mismatch. input array must have length
    equal to length(space) in its first dimension"
    K   = size(u, 2)

    # transform input
    sin = (ssp..., K)
    U   = _reshape(u, sin)

    # transform object
    ftr = plan_rfft(U, 1:D)

    # transform output
    V    = ftr * U
    sout = size(V)

    # output prototype
    M    = length(V) ÷ K
    sret = u isa AbstractMatrix ? (M, K) : (M,)
    v    = _reshape(V, sret)

    function fwd(v, u, p, t)
        U = _reshape(u, sin)
        V = _reshape(v, sout)
        mul!(V, ftr, U)

        v
    end

    function bwd(v, u, p, t)
        U = _reshape(u, sout)
        V = _reshape(v, sin)
        ldiv!(V, ftr, U)

        v
    end

    # form brfftplan, apply scaling. Look at rrule in AbstractFFTs
    function adj(v, u, p, t)
        U = _reshape(u, sout)
        V = _reshape(v, sin)

        v
    end

    ComplexT = if T isa Type{Float16}
        ComplexF16
    elseif T isa Type{Float32}
        ComplexF32
    else
        ComplexF64
    end

    ftransform = FunctionOperator(
                                  fwd;
                                  isinplace=true,
                                  T=ComplexT,
                                  size=(M,N),
    
                                  input_prototype=u,
                                  output_prototype=v,
    
                                  #op_adjoint= adj
                                  op_inverse = bwd
                                 )

    ftransform
end

## TODO - local system <-> global system
## global system for computation
## local system for plotting
# local_numbering(space::FourierSpace)
# global_numbering(space::FourierSpace)

###
# vector calculus
###

function massOp(space::FourierSpace{<:Any,1}, ::Galerkin)
    w = mass_matrix(space)
    DiagonalOperator(w)
end

###
# TODO - review gradOp(::FourierSpace) https://math.mit.edu/~stevenj/fft-deriv.pdf
# TODO   before writing vector calculus ops, transform operation on space
###

function gradOp(space::FourierSpace{<:Any,1})
    tr = transformOp(space)

    (k,) = modes(space)
    ik = im * DiagonalOperator(k)

    [
     tr \ ik * tr,
    ]
end

function hessianOp(space::FourierSpace{<:Any,1})
    tr = transformOp(space)

    (k,) = modes(space)
    ik2 = -DiagonalOperator(@. k * k)

    [
     tr \ ik2 * tr,
    ]
end

function laplaceOp(space::FourierSpace{<:Any,1}, ::Collocation)
    hessianOp(space) |> sum
end

function biharmonicOp(space::FourierSpace{<:Any,1})
    tr = transformOp(space)

    (k,) = modes(space)
    ik4 = DiagonalOperator(@. k * k)

    [
     tr \ ik4 * tr,
    ]
end

###
# interpolation operators
###

function interpOp(space1::FourierSpace{<:Any,1}, space2::FourierSpace{<:Any,1})
    tr1 = transformOp(space1)

    k1 = modes(space1)
    k2 = modes(space2)

    M = length.(k2) # output
    N = length.(k1) # input

    J = sparse(I, (M,N)) |> MatrixOperator

    tr1 \ J * tr1
end

transformOp(space::FourierSpace) = space.ftransform

###
# operators in transformed space
###

function gradOp(space::TransformedSpace{<:Any,1,<:FourierSpace})
    (k,) = modes(space)
    ik = DiagonalOperator(im * k)

    ik
end

function hessianOp(space::TransformedSpace{<:Any,1,<:FourierSpace})
    (k,) = modes(space)
    ik2 = DiagonalOperator(@. -k * k)

    ik2
end

#=
function advectionOp(vel::NTuple{D}, space::TransformedSpace{<:Any,D,<:FourierSpace}, discr::AbstractDiscretization) where{D}
    VV = [DiagonalOperator.(vel)...]

    tr = transformOp(space)

    VV_phys = tr \ VV

    MM = massOp(space, discr)
    DD = gradOp(space, discr)

    Dphys = tr \ DD

    Adv = _transp(VV_phys) * MM * Dphys


    -tr * Adv
end
=#

# interpolation
function interpOp(space1::TransformedSpace{<:Any,D,<:FourierSpace},
                  space2::TransformedSpace{<:Any,D,<:FourierSpace},
                 ) where{D}

    M = size(space2)[1] # output
    N = size(space1)[1] # input

    J = sparse(I, (M,N)) |> MatrixOperator

    J
end
#
