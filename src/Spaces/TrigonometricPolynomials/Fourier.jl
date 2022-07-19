#
# TODO - local system <-> global system
# global system for computation
# local system for plotting
# local_numbering(space::FourierSpace)
# global_numbering(space::FourierSpace)

"""
Fourier spectral space
"""
struct FourierSpace{
                    T,
                    D,
                    Tpts<:NTuple{D},
                    Tmds<:NTuple{D},
                    Tdom<:AbstractDomain{T,D},
                    Tgrid,
                    Tmodes,
                    Tmass,
                    Ttr,
                   } <: AbstractSpace{T,D}
    """ points """
    npoints::Tpts
    """ modes """
    nmodes::Tmds
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

function adapt_structure(to, space::FourierSpace)
    grid  = adapt_structure(to, space.grid)
    modes = adapt_structure(to, space.modes)
    mass_matrix = adapt_structure(to, space.mass_matrix)

    x = first(grid)
    T = eltype(x)

    npoints = space.npoints
    nmodes  = space.nmodes
    domain = T(space.domain)
    ftransform = form_transform(x, space)

    FourierSpace(
                 npoints, nmodes, domain, grid, modes,
                 mass_matrix, ftransform,
                )
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
    (L,) = lengths(domain)
    #""" reset deformation to map from [-π,π]^D """
    #ref_domain = reference_box(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO


    dz = L / n
    z  = range(start=-L/2, stop=L/2-dz, length=n) |> Array
    T  = eltype(z)

    FFTLIB = _fft_lib(z)
    k = FFTLIB.rfftfreq(n, 2π*n/L) |> Array

    npoints = (n,)
    nmodes  = (length(k),)
    domain = T(domain)
    grid = (z,)
    modes = (k,)
    mass_matrix = ones(T, n) * (2π/L)
    ftransform = nothing

    space = FourierSpace(
                         npoints, nmodes, domain, grid, modes,
                         mass_matrix, ftransform,
                        )

    space = make_transform(space, z)

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

function FourierSpace(nr::Integer, ns::Integer;
                      domain::AbstractDomain{<:Any,2}=FourierDomain(2),
                     )

    if !(domain isa BoxDomain)
        @error "Trigonometric polynomials work with logically rectangular domains"
    end

    domain = FourierDomain(2)
    (Lr, Ls) = lengths(domain)
    # reset deformation to map from [-π,π]^D
    #ref_domain = reference_box(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO

    dr = Lr / nr
    ds = Ls / ns
    zr = range(start=-Lr/2, stop=Lr/2-dr, length=nr) |> Array
    zs = range(start=-Lr/2, stop=Lr/2-ds, length=ns) |> Array

    FFTLIB = _fft_lib(zr)
    kr = FFTLIB.rfftfreq(nr, 2π*nr/Lr) |> Array
    ks = FFTLIB.fftfreq(ns, 2π*ns/Ls)  |> Array
    nkr = length(kr)
    nks = length(ks)

    r, s   = _vec.(ndgrid(zr, zs))
    kr, ks = _vec.(ndgrid(kr, ks))

    T  = eltype(r)

    npoints = (nr, ns)
    nmodes  = (nkr, nks)
    domain  = T(domain)
    grid    = (r, s)
    modes   = (kr, ks)
    mass_matrix = ones(T, nr * ns) * (2π/Lr) * (2π/Ls)
    ftransform  = nothing

    space = FourierSpace(
                         npoints, nmodes, domain, grid, modes,
                         mass_matrix, ftransform,
                        )

    space = make_transform(space, r)

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

###
# interface
###

Base.size(space::FourierSpace) = space.npoints
mode_size(space::FourierSpace) = space.nmodes
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
    K = size(u, 2)

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
        ldiv!(V, ftr, U) # TODO - check if fftlib caches inv(plan)

        v
    end

    ComplexT = if T isa Type{Float16}
        ComplexF16
    elseif T isa Type{Float32}
        ComplexF32
    else
        ComplexF64
    end

    FunctionOperator(
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
end

###
# operators in phys space
###

transformOp(space::FourierSpace) = space.ftransform

function truncationOp(space::FourierSpace, fracs=nothing)
    X = truncationOp(transform(space), fracs)

    if X isa IdentityOperator
        return IdentityOperator(space)
    end

    F = transformOp(space)

    F \ X * F
end

function massOp(space::FourierSpace, ::Galerkin)
    w = mass_matrix(space)
    DiagonalOperator(w)
end

function gradientOp(space::FourierSpace{<:Any,D}) where{D}
    sph = transform(space)  # transformed space
    DDh = gradientOp(sph)   # ∇ in transformed space

    F  = transformOp(space) # forward transform
    FF = [F for i=1:D]

    # https://github.com/vpuri3/PDEInterfaces.jl/issues/25
    FF .\ DDh .* FF # TODO - this is doing transform D times. should only be 1
end

function hessianOp(space::FourierSpace{<:Any,D}) where{D}
    sph  = transform(space)
    DD2h = hessianOp(sph)

    F  = transformOp(space)
    FF = [F for i=1:D]

    FF .\ DD2h .* FF
end

function laplaceOp(space::FourierSpace{<:Any,D}, discr::Collocation) where{D}
    sph = transform(space)
    D2h = laplaceOp(sph, discr)

    F = transformOp(space)

    F \ D2h * F
end

function biharmonicOp(space::FourierSpace{<:Any,D}) where{D}
    sph  = transform(space)
    DD4h = biharmonicOp(sph)

    F  = transformOp(space)
    FF = [F for i=1:D]

    FF .\ DD4h .* FF
end

function _fusedGradientTruncationOp(space::FourierSpace{<:Any,D},
                                    truncation_fracs=nothing,
                                   ) where{D}
    tspace = transform(space)

    F   = transformOp(space)
    Xh  = truncationOp(tspace, truncation_fracs)
    DDh = gradientOp(tspace)

    FF  = [F  for i=1:D]
    XXh = [Xh for i=1:D]

    FF .\ XXh .* DDh .* FF
end

function advectionOp(vels::NTuple{D},
                     space::FourierSpace{<:Any,D},
                     discr::AbstractDiscretization;
                     vel_update_funcs=nothing,
                     truncation_fracs=nothing,
                    ) where{D}

    VV = _pair_update_funcs(vels, vel_update_funcs)
    M  = massOp(space, discr)

    C = if M isa IdentityOperator # Collocation
        tspace = transform(space)

        F   = transformOp(space)
        Xh  = truncationOp(tspace, truncation_fracs)
        DDh = gradientOp(tspace)

        FF  = [F  for i=1:D]
        XXh = [Xh for i=1:D]

        VV' * (Diagonal(FF) \ Diagonal(XXh) * Diagonal(DDh)) * FF
    else # Galerkin
        X   = truncationOp(space, truncation_fracs)
        XDD = _fusedGradientTruncationOp(space, truncation_fracs)

        XX = Diagonal([X for i=1:D])
        MM = Diagonal([M for i=1:D])
        (XX*VV)' * MM * XDD
    end

    C
end

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

function truncationOp(space::TransformedSpace{<:Any,D,<:FourierSpace},
                      fracs::Union{NTuple{D,Number},Nothing}=nothing) where{D}

    fracs = fracs isa Nothing ? ([2//3 for d=1:D]) : fracs

    if isone(prod(fracs))
        return IdentityOperator(space)
    end

    ns = size(space)

    a = ones(Bool, ns)
    for d=1:D
        n = ns[d]
        frac = fracs[d]

        idx = if d == 1
            cut = (n-1)*frac |> round |> Int

            n-cut+1 : n
        else
            mid = n/2 + 1 |> round |> Int
            cut = (n-1)/2 * frac |> round |> Int

            mid-cut+1 : mid+cut-1
        end

        a[(Colon() for i=1:d-1)..., idx, (Colon() for i=d+1:D)...] .= false
    end
    a = points(space)[1] isa CUDA.CuArray ? gpu(a) : a

    DiagonalOperator(_vec(a))
end

function gradientOp(space::TransformedSpace{<:Any,D,<:FourierSpace}) where{D}
    ks = points(space)
    ns = size(transform(space))

    # https://math.mit.edu/~stevenj/fft-deriv.pdf
    iks = [@. im*ks[i] for i=1:D]
    for i=1:D iseven(ns[i]) && CUDA.@allowscalar iks[i][end] = 0  end

    DiagonalOperator.(iks)
end

function hessianOp(space::TransformedSpace{<:Any,D,<:FourierSpace}) where{D}
    ks = points(space)
    ik2s = [@. -ks[i]^2 for i=1:D]

    DiagonalOperator.(ik2s)
end

function laplaceOp(space::TransformedSpace{<:Any,D,FourierSpace}, ::Collocation) where{D}
    ks = points(space)
    ik2 = [@. -ks[i]^2 for i=1:D] |> sum

    DiagonalOperator.(ik2)
end

function biharmonicOp(space::TransformedSpace{<:Any,D,<:FourierSpace}) where{D}
    ks = points(space)
    ik4s = [@. ks[i]^4 for i=1:D]

    DiagonalOperator.(ik4s)
end

function advectionOp(vels::NTuple{D},
                     tspace::TransformedSpace{<:Any,D,<:FourierSpace},
                     discr::AbstractDiscretization;
                     vel_update_funcs=nothing,
                     truncation_fracs=nothing,
                    ) where{D}

    @error "this method has problems. fix later"
    VV = _pair_update_funcs(vels, vel_update_funcs)

    DDh = gradientOp(tspace, discr)

    # physical space
    space = transform(tspace)

    F  = transformOp(space)
    M  = massOp(space, discr)
    Xh = truncationOp(tspace, truncation_fracs)

#   C = if M isa IdentityOperator
#   else
#   end

    MM  = Diagonal([M  for i=1:D])
    FF  = Diagonal([F  for i=1:D])
    XXh = Diagonal([Xh for i=1:D])

    VV_phys = FF \ XXh * VV
    DD_phys = FF \ XXh * DDh

    # problem: V * (Dh * û0) == 0 <-- problem

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
