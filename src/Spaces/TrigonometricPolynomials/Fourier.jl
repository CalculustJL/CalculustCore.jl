#
"""
with stuff in physical space. operators use FFT
"""
struct FourierSpace{
                    T,
                    D,
                    Tdom<:AbstractDomain{T,D},
                    Tgrid,
                    Tmodes,
                    Tmass,
                    Ttr,
                   } <: AbstractSpace{T,D}
    """ Domain """
    domain::Tdom
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


    dx = L / n
    x  = range(start=-L/2, stop=L/2-dx, length=n) |> Array
    T  = eltype(x)

    FFTLIB = FFTW #_fft_lib(x)
    k = FFTLIB.rfftfreq(n, 2π*n/L) |> Array

    domain = T(domain)
    grid = (x,)
    modes = (k,)
    mass_matrix = ones(T, n) * (2π/L)
    ftransform = nothing

    space = FourierSpace(
                         domain, grid, modes,
                         mass_matrix, ftransform,
                        )

    space = make_transform(space, x)

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

# TODO - just @functor FourierSpace instead
function Adapt.adapt_storage(::LuxCUDAAdaptor, space::FourierSpace)
    grid = CUDA.cu(space.grid)
    modes = CUDA.cu(space.modes)
    mass_matrix = CUDA.cu(space.mass_matrix)

    x = first(grid)
    T = eltype(x)

    domain = T(space.domain)
    ftransform = form_transform(x, space)

    FourierSpace(
                 domain, grid, modes,
                 mass_matrix, ftransform,
                )
end

###
# interface
###

Base.size(space::FourierSpace) = points(space) |> first |> size
domain(space::FourierSpace) = space.domain
points(space::FourierSpace) = space.grid
function quadratures(space::FourierSpace{<:Any,1})
    x = points(space) |> first
    w = mass_matrix(space)

    ((x, w),)
end
mass_matrix(space::FourierSpace) = space.mass_matrix
modes(space::FourierSpace) = space.modes

function form_transform(u::AbstractVecOrMat{T}, space::FourierSpace{<:Any,D};
                        p=nothing, t=zero(T)) where{T,D}

    ssp = size(space)
    N   = length(space)

    @assert size(u, 1) == N "size mismatch. input array must have length
    equal to length(space) in its first dimension"
    K   = size(u, 2)

    # transform input shape
    sin = (ssp..., K)
    U   = _reshape(u, sin)

    # transform object
    FFTLIB = _fft_lib(u)
    ftr = FFTLIB.plan_rfft(U, 1:D)

    # transform output shape
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

                                  op_inverse=bwd,
                                  op_adjoint=bwd,
                                  op_adjoint_inverse=fwd,

                                  p=p,
                                  t=t,
                                 )

    ftransform
end

###
# Operators
###

transformOp(space::FourierSpace) = space.ftransform

function truncationOp(space::FourierSpace{<:Any,1}, frac=nothing)
    X = truncationOp(transform(space), frac)

    if X isa IdentityOperator
        return IdentityOperator(space)
    end

    F = transformOp(space)

    F \ X * F
end

function truncationOp(space::TransformedSpace{<:Any,1,<:FourierSpace}, frac=nothing)

    frac = frac isa Nothing ? 2//3 : frac
    if isone(frac)
        return IdentityOperator(space)
    end

    (n,) = length.(points(space))

    a = begin
        a = [true for i=1:n]
        m = n * frac |> round |> Int
        a[m:n] .= false

        points(space)[1] isa CUDA.CuArray ? gpu(a) : a
    end

    DiagonalOperator(a)
end

###
# vector calculus
###

## TODO - local system <-> global system
## global system for computation
## local system for plotting
# local_numbering(space::FourierSpace)
# global_numbering(space::FourierSpace)

###
# TODO - FourierSpace operators - review https://math.mit.edu/~stevenj/fft-deriv.pdf
###

function massOp(space::FourierSpace{<:Any,1}, ::Galerkin)
    w = mass_matrix(space)
    DiagonalOperator(w)
end

function gradientOp(space::FourierSpace{<:Any,D}) where{D}
    sph = transform(space)  # transformed space
    DDh = gradientOp(sph)   # ∇ in transformed space

    F  = transformOp(space) # forward transform
    FF = [F for i=1:D]

    FF .\ DDh .* FF
end

function hessianOp(space::FourierSpace{<:Any,D}) where{D}
    sph  = transform(space)
    DD2h = hessianOp(sph)

    F  = transformOp(space)
    FF = [F for i=1:D]

    FF .\ DD2h .* FF
end

function biharmonicOp(space::FourierSpace{<:Any,D}) where{D}
    sph  = transform(space)
    DD4h = biharmonicOp(sph)

    F  = transformOp(space)
    FF = [F for i=1:D]

    FF .\ DD4h .* FF
end

function _fusedGradientTruncationOp(space::FourierSpace{<:Any,D},
                                    truncation_frac=nothing,
                                   ) where{D}
    tspace = transform(space)

    F   = transformOp(space)
    Xh  = truncationOp(tspace, truncation_frac)
    DDh = gradientOp(tspace)

    FF  = [F  for i=1:D]
    XXh = [Xh for i=1:D]

    FF .\ XXh .* DDh .* FF
end

function advectionOp(vels::NTuple{D},
                     space::FourierSpace{<:Any,D},
                     discr::AbstractDiscretization;
                     vel_update_funcs=nothing,
                     truncation_frac=nothing,
                    ) where{D}

    VV = _pair_update_funcs(vels, vel_update_funcs)
    M  = massOp(space, discr)

    C = if M isa IdentityOperator # Collocation
        tspace = transform(space)

        F   = transformOp(space)
        Xh  = truncationOp(tspace, truncation_frac)
        DDh = gradientOp(tspace)

        FF  = [F  for i=1:D]
        XXh = [Xh for i=1:D]

        VV' * (Diagonal(FF) \ Diagonal(XXh) * Diagonal(DDh)) * FF
    else # Galerkin
        X   = truncationOp(space, truncation_frac)
        XDD = _fusedGradientTruncationOp(space, truncation_frac)

        XX = Diagonal([X for i=1:D])
        MM = Diagonal([M for i=1:D])
        (XX*VV)' * MM * XDD
    end

    C
end

###
# interpolation operators
###

function interpOp(space1::FourierSpace{<:Any,1}, space2::FourierSpace{<:Any,1})
    F1   = transformOp(space1)
    F2   = transformOp(space2)
    sp1h = transform(space1)
    sp2h = transform(space2)

    J = interpOp(sp1h, sp2h)

    F2 \ J * F1
end

###
# operators in transformed space
###

function gradientOp(space::TransformedSpace{<:Any,1,<:FourierSpace})
    (k,) = points(space)
    ik = DiagonalOperator(im * k)

    [
     ik,
    ]
end

function hessianOp(space::TransformedSpace{<:Any,1,<:FourierSpace})
    (k,) = points(space)
    ik2 = DiagonalOperator(@. -k * k)

    [
     ik2,
    ]
end

function biharmonicOp(space::TransformedSpace{<:Any,1,FourierSpace})
    (k,) = points(space)
    ik4 = DiagonalOperator(@. k^4)

    [
     ik4,
    ]
end

function advectionOp(vels::NTuple{D},
                     tspace::TransformedSpace{<:Any,D,<:FourierSpace},
                     discr::AbstractDiscretization;
                     vel_update_funcs=nothing,
                     truncation_frac=nothing,
                    ) where{D}

    @error "this method has problems. fix later"
    VV = _pair_update_funcs(vels, vel_update_funcs)

    DDh = gradientOp(tspace, discr)

    # physical space
    space = transform(tspace)

    F  = transformOp(space)
    M  = massOp(space, discr)
    Xh = truncationOp(tspace, truncation_frac)

#   C = if M isa IdentityOperator
#   else
#   end

    MM  = Diagonal([M  for i=1:D])
    FF  = Diagonal([F  for i=1:D])
    XXh = Diagonal([Xh for i=1:D])

    VV_phys = FF \ XXh * VV
    DD_phys = FF \ XXh * DDh

    # V * (Dh * û0) == 0 <-- problem

    VV_phys' * MM * DD_phys # <- this is zero
#   transpose(VV) * XXh * FF * MM * DD_phys
end

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
