
###
# Lagrange polynomial function spaces
###

""" Lagrange polynomial spectral space """
struct LagrangePolynomialSpace{T,
                               D,
                               Tpts,
                               Tdom<:AbstractDomain{T,D},
                               Tquad,
                               Tgrid,
                               Tmass,
                               Tderiv,
                               Tloc,
#                              Tglo
                              } <: AbstractSpectralSpace{T,D}
    """ Domain """
    dom::Tdom
    """ size tuple """
    npoints::Tpts
    """ quadratures """
    quads::Tquad
    """ grid points """
    points::Tgrid
    """ mass matrix """
    mass_matrix::Tmass
    """ derivative matrices """
    deriv_mats::Tderiv
    """ local numbering """
    loc_num::Tloc
    """ global numbering """
#   glo_num::Tglo
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

    dom = T(domain)
    npoints = (n,)
    quads = ((z, w),)
    points = _vec.((z,))
    mass_matrix = _vec(w)
    deriv_mats = (D,)
    local_numbering = _reshape(1:prod(npoints), npoints)

    space = LagrangePolynomialSpace(
                                    domain, npoints, quads, points,
                                    mass_matrix, deriv_mats, 
                                    local_numbering,
                                   )

    domain isa DeformedDomain ? deform(space, mapping) : space
end

GaussLobattoLegendre1D(args...; kwargs...) =
    LagrangePolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendre1D(args...; kwargs...) =
    LagrangePolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebychev1D(args...; kwargs...) =
    LagrangePolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)

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
    quads = ((zr, wr), (zs, ws),)
    points = _vec.((r, s,))
    mass_matrix = _vec(wr * ws')
    deriv_mats = (Dr, Ds,)
    local_numbering = _reshape(1:prod(npoints), npoints)

    space = LagrangePolynomialSpace(
                                    domain, npoints, quads, points,
                                    mass_matrix, deriv_mats,
                                    local_numbering,
                                   )

    domain isa DeformedDomain ? deform(space, mapping) : space
end

GaussLobattoLegendre2D(args...; kwargs...) =
    LagrangePolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendre2D(args...; kwargs...) =
    LagrangePolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebychev2D(args...; kwargs...) =
    LagrangePolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)

### abstract interface

Base.size(space::LagrangePolynomialSpace) = space.npoints

function Plots.plot(u, space::LagrangePolynomialSpace{<:Number,1})

    (x,) = grid(space)
    plt  = plot(x, u, legend=false)
end

function Plots.plot(u, space::LagrangePolynomialSpace{<:Number,2}; a=45, b=60)

    npts = size(space)
    (x,y,) = grid(space)

    u = _reshape(u, npts)
    x = _reshape(x, npts)
    y = _reshape(y, npts)

    plt = plot(x, y, u, legend=false, c=:grays, camera=(a,b))
    plt = plot!(x', y', u', legend=false, c=:grays, camera=(a,b))

    plt
end

grid(space::LagrangePolynomialSpace) = space.points
domain(space::LagrangePolynomialSpace) = space.dom
quadratures(space::LagrangePolynomialSpace) = space.quads
local_numbering(space::LagrangePolynomialSpace) = space.loc_num

function global_numbering(space::AbstractSpace)
    dom = domain(space)
    loc_num = local_numbering(space)
end

function boundary_nodes(space::LagrangePolynomialSpace)
    D = dims(space)
    npoints = size(space)
    loc_num = local_numbering(space)

    indices = ()
    for i=1:D
        n = npoints[i]
        range_lower = ([1:npoints[j] for j=1:i-1]..., 1, [1:npoints[j] for j=i+1:D]...)
        range_upper = ([1:npoints[j] for j=1:i-1]..., n, [1:npoints[j] for j=i+1:D]...)
        indices = (indices..., loc_num[range_lower...])
        indices = (indices..., loc_num[range_upper...])
    end

    indices
end

### vector calculus ops

function massOp(space::LagrangePolynomialSpace)
    @unpack mass_matrix = space

    DiagonalOperator(mass_matrix)
end

function gradOp(space::LagrangePolynomialSpace{<:Number,1})
    (Dr,) = space.deriv_mats

    Dx = MatrixOperator(Dr)

    DD = AbstractSciMLOperator[Dx]
end

function gradOp(space::LagrangePolynomialSpace{<:Number,2})
    (nr, ns) = space.npoints
    (Dr, Ds) = space.deriv_mats

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()

    Dx = ⊗(Dr, Is)
    Dy = ⊗(Ir, Ds)

    DD = AbstractSciMLOperator[Dx
                               Dy]
end

function gradOp(space::LagrangePolynomialSpace{<:Number,3})
    (Dr, Ds, Dt) = space.deriv_mats
    (nr, ns, nt) = space.npoints

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()
    It = IdentityOperator{nt}()

    Dx = ⊗(Dr, Is, It)
    Dy = ⊗(Ir, Ds, It)
    Dz = ⊗(Ir, Is, Dt)

    DD = AbstractSciMLOperator[Dx
                               Dy
                               Dz]
end

### interpolation operators

function interpOp(space1::LagrangePolynomialSpace{<:Number,1},
                  space2::LagrangePolynomialSpace{<:Number,1})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    J = lagrange_interp_mat(r2, r1) # from 1 to 2

    MatrixOperator(J)
end

function interpOp(space1::LagrangePolynomialSpace{<:Number,2},
                  space2::LagrangePolynomialSpace{<:Number,2})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    s1, _ = space1.quads[2]
    s2, _ = space2.quads[2]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)

    ⊗(Jr, Js)
end

function interpOp(space1::LagrangePolynomialSpace{<:Number,3},
                  space2::LagrangePolynomialSpace{<:Number,3})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    s1, _ = space1.quads[2]
    s2, _ = space2.quads[2]

    t1, _ = space1.quads[3]
    t2, _ = space2.quads[3]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)
    Jt = lagrange_interp_mat(t2, t1)

    ⊗(Jr, Js, Jt)
end
#
