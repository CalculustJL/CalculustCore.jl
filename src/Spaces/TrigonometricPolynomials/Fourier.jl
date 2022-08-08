#
"""
Fourier spectral space
"""
struct FourierSpace{
                    T,
                    D,
                    Tnpoints<:NTuple{D,Integer},
                    Tnfreqs<:NTuple{D,Integer},
                    Tdom<:AbstractDomain{T,D},
                    Tgrid<:NTuple{D,AbstractArray{T}},
                    Tfreqs<:NTuple{D,AbstractArray{T}},
                    Tmass_mat<:AbstractArray{T},
                    Tftransform,
                   } <: AbstractSpace{T,D}
    """ # points """
    npoints::Tnpoints
    """ # freq size """
    nfreqs::Tnfreqs
    """ Domain """
    dom::Tdom
    """ grid points """
    grid::Tgrid
    """ frequencies """
    freqs::Tfreqs
    """ mass matrix """
    mass_mat::Tmass_mat
    """ forward transform `mul!(û, T , u)` """
    ftransform::Tftransform
end

function (::Type{T})(space::FourierSpace) where{T<:Number}
    npoints = size(space)
    nfreqs  = mode_size(space)
    dom     = T(domain(space))

    grid  = Tuple(T.(x) for x in points(space))
    freqs = Tuple(T.(k) for k in modes(space))

    mass_mat = T.(mass_matrix(space))
    ftransform = transformOp(space)

    space = FourierSpace(
                         npoints, nfreqs, dom, grid, freqs,
                         mass_mat, nothing,
                        )

    p   = nothing # TODO
    u0  = T.(ftransform.cache[1])

    make_transform(space, u0; p=p)
end

function adapt_structure(to, space::FourierSpace)
    grid  = adapt_structure(to, points(space))
    freqs = adapt_structure(to, modes(space))
    mass_mat = adapt_structure(to, mass_matrix(space))

    x = first(grid)
    T = eltype(x)

    npoints = adapt_structure(to, size(space))
    nfreqs  = adapt_structure(to, mode_size(space))
    dom     = adapt_structure(to, T(domain(space)))
    ftransform = transformOp(space)

    space = FourierSpace(
                         npoints, nfreqs, dom, grid, freqs,
                         mass_mat, nothing,
                        )

    p   = adapt_structure(to, ftransform.p)
    u0  = adapt_structure(to, ftransform.cache[1])

    make_transform(space, u0; p=p)
end

function FourierSpace(n::Integer;
                      domain::AbstractDomain{<:Any,1}=FourierDomain(1),
                     )

    dom = if domain isa IntervalDomain
        BoxDomain(domain)
    elseif domain isa BoxDomain
        domain
    else
        @error "Trigonometric polynomials work with logically rectangular domains"
    end

    dom = FourierDomain(1)
    (L,) = lengths(dom)
    #""" reset deformation to map from [-π,π]^D """
    #ref_dom = reference_box(2)
    #dom = ref_dom # map_from_ref(dom, ref_dom) # TODO


    dz = L / n
    z  = range(start=-L/2, stop=L/2-dz, length=n) |> Array
    T  = eltype(z)

    FFTLIB = _fft_lib(z)
    k = FFTLIB.rfftfreq(n, 2π*n/L) |> Array

    npoints = (n,)
    nfreqs  = (length(k),)
    dom = T(dom)
    grid = (z,)
    freqs = (k,)
    mass_mat = ones(T, n) * (2π/L)
    ftransform = nothing

    space = FourierSpace(
                         npoints, nfreqs, dom, grid, freqs,
                         mass_mat, ftransform,
                        )

    space = make_transform(space, z)

    dom isa Domains.DeformedDomain ? deform(space, mapping) : space
end

function FourierSpace(nr::Integer, ns::Integer;
                      domain::AbstractDomain{<:Any,2}=FourierDomain(2),
                     )

    dom = if domain isa BoxDomain
        domain
    else
        @error "Trigonometric polynomials work with logically rectangular domains"
    end

    dom = FourierDomain(2)
    (Lr, Ls) = lengths(dom)
    # reset deformation to map from [-π,π]^D
    #ref_dom = reference_box(2)
    #dom = ref_dom # map_from_ref(dom, ref_dom) # TODO

    dr = Lr / nr
    ds = Ls / ns
    zr = range(start=-Lr/2, stop=Lr/2-dr, length=nr) |> Array
    zs = range(start=-Lr/2, stop=Lr/2-ds, length=ns) |> Array

    FFTLIB = _fft_lib(zr)
    kr = FFTLIB.rfftfreq(nr, 2π*nr/Lr) |> Array
    ks = FFTLIB.fftfreq(ns, 2π*ns/Ls)  |> Array
    nkr = length(kr)
    nks = length(ks)

    r, s   = vec.(ndgrid(zr, zs))
    kr, ks = vec.(ndgrid(kr, ks))

    T  = eltype(r)

    npoints = (nr, ns)
    nfreqs  = (nkr, nks)
    dom  = T(dom)
    grid    = (r, s)
    freqs   = (kr, ks)
    mass_mat = ones(T, nr * ns) * (2π/Lr) * (2π/Ls)
    ftransform  = nothing

    space = FourierSpace(
                         npoints, nfreqs, dom, grid, freqs,
                         mass_mat, ftransform,
                        )

    space = make_transform(space, r)

    dom isa Domains.DeformedDomain ? deform(space, mapping) : space
end

###
# interface
###

Base.size(space::FourierSpace) = space.npoints
mode_size(space::FourierSpace) = space.nfreqs
domain(space::FourierSpace) = space.dom
points(space::FourierSpace) = space.grid
function quadratures(space::FourierSpace{<:Any,1})
    x = points(space) |> first
    w = mass_matrix(space)

    ((x, w),)
end
mass_matrix(space::FourierSpace) = space.mass_mat
modes(space::FourierSpace) = space.freqs

function form_transform(
                        space::FourierSpace{<:Any,D},
                        u::Union{Nothing,AbstractVecOrMat}=nothing;
                        p=nothing,
                        t::Union{Real,Nothing}=nothing,
                       ) where{D}

    # sinput : (N,K) - input size
    # soutput: (M,K) - output size
    #
    # sin : (n1,...,nd) - input size to FFT st n1*...*nd=N
    # sout: (m1,...,md) - output size to FFT st m1*...*md=M

    u = u isa Nothing ? points(space) |> first : u
    T = eltype(u)
    t = zero(T)

    sinput = size(u)
    sspace = size(space)
    N   = length(space)

    @assert size(u, 1) == N "size mismatch. input array must have length
    $(length(space)) in its first dimension"
    K = size(u, 2)

    # transform input shape
    sin = (sspace..., K)
    U   = reshape(u, sin)

    # transform object
    FFTLIB = _fft_lib(u)
    ftr = FFTLIB.plan_rfft(U, 1:D)

    # transform output shape
    V    = ftr * U
    sout = size(V)

    # output prototype
    M = length(V) ÷ K
    soutput = u isa AbstractMatrix ? (M, K) : (M,)
    v    = reshape(V, soutput)

    # in-place
    function fwd(v, u, p, t)
        U = reshape(u, sin)
        V = reshape(v, sout)
        mul!(V, ftr, U)

        v
    end

    function bwd(v, u, p, t)
        U = reshape(u, sout)
        V = reshape(v, sin)
        ldiv!(V, ftr, U)

        v
    end

    # out-of-place
    function fwd(u, p, t)
        U = reshape(u, sin)
        V = ftr * U

        reshape(V, soutput)
    end

    function bwd(u, p, t)
        U = reshape(u, sout)
        V = ftr \ U

        reshape(V, sinput)
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
                     outofplace=true,
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

function interpOp(space1::FourierSpace{T1,1}, space2::FourierSpace{T2,1}) where{T1,T2}
    F1   = transformOp(space1)
    F2   = transformOp(space2)
    sp1h = transform(space1)
    sp2h = transform(space2)

    J = interpOp(sp1h, sp2h)

    N1 = length(space1)
    N2 = length(space2)
    λ = T1(N2 / N1) # TODO - verify in AbstractFFTs documentation

    λ * F1 \ J * F2
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
            cut = round(Int, n*frac, RoundUp)

            cut : n
        else
            cut = round(Int, n ÷ 2 * frac, RoundUp)

            cut : n-cut
        end

        a[(Colon() for i=1:d-1)..., idx, (Colon() for i=d+1:D)...] .= false
    end
    a = points(space)[1] isa CUDA.CuArray ? gpu(a) : a

    DiagonalOperator(vec(a))
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
                     tspace::TransformedSpace{T,D,<:FourierSpace},
                     discr::AbstractDiscretization;
                     vel_update_funcs=nothing,
                     truncation_fracs=nothing,
                    ) where{T,D}

    space = transform(tspace)

    F  = transformOp(space)
    Xh = truncationOp(tspace, truncation_fracs)
    M  = massOp(space, discr)

    VV = begin
        trunc = cache_operator(F \ Xh, vels[1])
        vel_funcs = ComposedUpdateFunction.((trunc,), vel_update_funcs, deepcopy(vels))

        _pair_update_funcs((F,) .\ vels, vel_funcs)
    end

    MM  = Diagonal([M  for i=1:D])
    XXh = Diagonal([Xh for i=1:D])
    FFi = Diagonal([inv(F)  for i=1:D])
    DDh = gradientOp(tspace, discr)

    F * (VV' * MM * FFi * XXh * DDh)
end

# interpolation
function interpOp(space1::TransformedSpace{<:Any,D,<:FourierSpace},
                  space2::TransformedSpace{<:Any,D,<:FourierSpace},
                 ) where{D}

    Ms = size(space1) # output
    Ns = size(space2) # input

    sz = Tuple(zip(Ms,Ns))
    Js = Tuple(sparse(I, sz[i]) for i=1:D)

    ⊗(Js...)
end
#
