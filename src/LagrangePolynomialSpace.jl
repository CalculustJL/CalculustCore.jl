#
###
# Lagrange polynomial function spaces
###

""" Lagrange polynomial spectral space """
struct LagrangePolynomialSpace{T,D,Tpts,
                               Tdom<:AbstractDomain{T,D},
                               Tquad,Tgrid,
                               Tmass<:AbstractField{T,D},
                               Tderiv,
#                              Tl2g,
                              } <: AbstractSpectralSpace{T,D}
    domain::Tdom
    npoints::Tpts
    quadratures::Tquad
    grid::Tgrid
    mass_matrix::Tmass
    deriv_mats::Tderiv
#   local2global::Tl2g # TODO
end

function LagrangePolynomialSpace(n::Integer;
        domain::AbstractDomain{<:Number,1}=reference_box(1),
        quadrature = gausslobatto,
        T = Float64,
       )

    if domain isa IntervalDomain
        domain = BoxDomain(domain)
    elseif !(domain isa BoxDomain)
        @error "spectral polynomials work with logically rectangular domains"
    end

    #""" reset deformation to map from [-1,1]^D """
    #ref_domain = reference_box(1)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO
    ## change domain eltype

    z, w = quadrature(n)

    z = T.(z)
    w = T.(w)

    D = lagrange_deriv_mat(z)

    domain = T(domain)
    npoints = (n,)
    quadratures = ((z, w),)
    grid = (Field(z),)
    mass_matrix = Field(w)
    deriv_mats = (D,)

    space = LagrangePolynomialSpace(
                                    domain, npoints, quadratures, grid,
                                    mass_matrix, deriv_mats, 
                                   )

    domain isa DeformedDomain ? deform(space, mapping) : space
end

GaussLobattoLegendre1D(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendre1D(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebychev1D(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)

function LagrangePolynomialSpace(nr::Integer, ns::Integer;
        domain::AbstractDomain{<:Number,2}=reference_box(2),
        quadrature = gausslobatto,
        T = Float64,
       )

    if !(domain isa BoxDomain)
        @error "spectral polynomials work with logically rectangular domains"
    end

    #""" reset deformation to map from [-1,1]^D """
    #ref_domain = reference_box(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO

    zr, wr = quadrature(nr)
    zs, ws = quadrature(ns)

    zr, wr = T.(zr), T.(wr)
    zs, ws = T.(zs), T.(ws)

    r, s = ndgrid(zr,zs)

    Dr = lagrange_deriv_mat(zr)
    Ds = lagrange_deriv_mat(zs)

    domain = T(domain)
    npoints = (nr, ns,)
    quadratures = ((zr, wr), (zs, ws),)
    grid = (Field(r), Field(s))
    mass_matrix = Field(wr * ws')
    deriv_mats = (Dr, Ds,)

    space = LagrangePolynomialSpace(
                                    domain, npoints, quadratures, grid,
                                    mass_matrix, deriv_mats, 
                                   )

    domain isa DeformedDomain ? deform(space, mapping) : space
end

GaussLobattoLegendre2D(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendre2D(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebychev2D(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)

### abstract interface

get_grid(space::LagrangePolynomialSpace) = space.grid
get_domain(space::LagrangePolynomialSpace) = space.domain
numpoints(space::LagrangePolynomialSpace) = space.npoints

function boundary_nodes(space::LagrangePolynomialSpace{<:Number,D}) where{D}
    npoints = numpoints(space)
    indices = []
    for i=1:D
        n = npoints[i]
        range_lower = ([Colon() for j=1:i-1]..., 1, [Colon() for j=i+1:D]...)
        range_upper = ([Colon() for j=1:i-1]..., n, [Colon() for j=i+1:D]...)
        push!(indices, range_lower)
        push!(indices, range_upper)
    end

    indices
end

### vector calculus ops

function massOp(space::LagrangePolynomialSpace)
    @unpack mass_matrix = space

    DiagonalOp(mass_matrix)
end

function gradOp(space::LagrangePolynomialSpace{<:Number,1})
    (Dr,) = space.deriv_mats

    Dx = MatrixOp(Dr)

    DD = AbstractOperator{<:Number,1}[Dx]
    reshape(DD, 1, 1)
end

function gradOp(space::LagrangePolynomialSpace{<:Number,2})
    (nr, ns) = space.npoints
    (Dr, Ds) = space.deriv_mats

    Ir = I #sparse(I, nr, nr)
    Is = I #sparse(I, ns, ns)

    Dx = TensorProductOp2D(Dr, Is)
    Dy = TensorProductOp2D(Ir, Ds)

    DD = AbstractOperator{<:Number,2}[Dx
                                      Dy]
    reshape(DD, 2, 1)
end

function gradOp(space::LagrangePolynomialSpace{<:Number,3})
    (Dr, Ds, Dt) = space.deriv_mats
    (nr, ns, nt) = space.npoints

    Ir = I #sparse(I, nr, nr)
    Is = I #sparse(I, ns, ns)
    It = I #sparse(I, nt, nt)

    Dx = TensorProductOp3D(Dr, Is, It)
    Dy = TensorProductOp3D(Ir, Ds, It)
    Dz = TensorProductOp3D(Ir, Is, Dt)

    DD = AbstractOperator{<:Number,2}[Dx
                                      Dy
                                      Dz]
    reshape(DD, 3, 1)
end

### interpolation operators

function interpOp(space1::LagrangePolynomialSpace{<:Number,1},
                  space2::LagrangePolynomialSpace{<:Number,1})
    r1, _ = space1.quadratures[1]
    r2, _ = space2.quadratures[1]

    J = lagrange_interp_mat(r2, r1) # from 1 to 2

    MatrixOp(J)
end

function interpOp(space1::LagrangePolynomialSpace{<:Number,2},
                  space2::LagrangePolynomialSpace{<:Number,2})
    r1, _ = space1.quadratures[1]
    r2, _ = space2.quadratures[1]

    s1, _ = space1.quadratures[2]
    s2, _ = space2.quadratures[2]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)

    TensorProductOp2D(Jr, Js)
end

function interpOp(space1::LagrangePolynomialSpace{<:Number,3},
                  space2::LagrangePolynomialSpace{<:Number,3})
    r1, _ = space1.quadratures[1]
    r2, _ = space2.quadratures[1]

    s1, _ = space1.quadratures[2]
    s2, _ = space2.quadratures[2]

    t1, _ = space1.quadratures[3]
    t2, _ = space2.quadratures[3]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)
    Jt = lagrange_interp_mat(t2, t1)

    TensorProductOp3D(Jr, Js, Jt)
end
#
