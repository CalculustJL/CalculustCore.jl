#
include("LagrangeMatrices.jl")

###
# Lagrange polynomial function spaces
###

"""
Lagrange polynomial spectral space
"""
struct LagrangePolynomialSpace{T,
                               D,
                               Tpts<:NTuple{D},
                               Tdom<:AbstractDomain{T,D},
                               Tquad,
                               Tgrid,
                               Tmass,
                               Tderiv,
                               Tloc,
#                              Tglo
                              } <: AbstractSpace{T,D}
    """ size """
    npoints::Tpts
    """ Domain """
    domain::Tdom
    """ quadratures """
    quads::Tquad
    """ grid points """
    grid::Tgrid
    """ mass matrix """
    mass_matrix::Tmass
    """ derivative matrices """
    deriv_mats::Tderiv
    """ local numbering """
    loc_num::Tloc
#   """ global numbering """
#   glo_num::Tglo
end

function LagrangePolynomialSpace(n::Integer;
        domain::AbstractDomain{<:Any,1}=ChebychevDomain(1),
        quadrature = gausslobatto,
        T = Float64,
       )

    if domain isa IntervalDomain
        domain = BoxDomain(domain)
    elseif !(domain isa BoxDomain)
        @error "spectral polynomials work with logically rectangular domains"
    end

    #""" reset deformation to map from [-1,1]^D """
    #ref_domain = ChebychevDomain(1)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO
    ## change domain eltype

    z, w = quadrature(n)

    z = T.(z)
    w = T.(w)

    D = lagrange_deriv_mat(z)

    npoints = (n,)
    domain = T(domain)
    quads  = ((z, w),)
    grid   = _vec.((z,))
    mass_matrix = _vec(w)
    deriv_mats = (D,)
    local_numbering = _reshape(1:prod(npoints), npoints)

    space = LagrangePolynomialSpace(
                                    npoints, domain, quads, grid,
                                    mass_matrix, deriv_mats, 
                                    local_numbering,
                                   )

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

function LagrangePolynomialSpace(nr::Integer, ns::Integer;
        domain::AbstractDomain{<:Number,2}=ChebychevDomain(2),
        quadrature = gausslobatto,
        T = Float64,
       )

    if !(domain isa BoxDomain)
        @error "spectral polynomials work with logically rectangular domains"
    end

    #""" reset deformation to map from [-1,1]^D """
    #ref_domain = ChebychevDomain(2)
    #domain = ref_domain # map_from_ref(domain, ref_domain) # TODO

    zr, wr = quadrature(nr)
    zs, ws = quadrature(ns)

    zr, wr = T.(zr), T.(wr)
    zs, ws = T.(zs), T.(ws)

    r, s = ndgrid(zr,zs)

    Dr = lagrange_deriv_mat(zr)
    Ds = lagrange_deriv_mat(zs)

    npoints = (nr, ns,)
    domain = T(domain)
    quads = ((zr, wr), (zs, ws),)
    grid = _vec.((r, s,))
    mass_matrix = _vec(wr * ws')
    deriv_mats = (Dr, Ds,)
    local_numbering = _reshape(1:prod(npoints), npoints)

    space = LagrangePolynomialSpace(
                                    npoints, domain, quads, grid,
                                    mass_matrix, deriv_mats,
                                    local_numbering,
                                   )

    domain isa Domains.DeformedDomain ? deform(space, mapping) : space
end

GaussLobattoLegendreSpace(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendreSpace(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebychevSpace(args...; kwargs...) = LagrangePolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)

###
# interface
###

Base.size(space::LagrangePolynomialSpace) = space.npoints

domain(space::LagrangePolynomialSpace) = space.domain
points(space::LagrangePolynomialSpace) = space.grid
quadratures(space::LagrangePolynomialSpace) = space.quads
mass_matrix(space::LagrangePolynomialSpace) = DiagonalOperator(space.mass_matrix)
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

###
# vector calculus ops
###

function massOp(space::LagrangePolynomialSpace, ::Galerkin)
    @unpack mass_matrix = space

    DiagonalOperator(mass_matrix)
end

function gradientOp(space::LagrangePolynomialSpace{<:Number,1})
    (Dr,) = space.deriv_mats

    Dx = MatrixOperator(Dr)

    DD = [Dx,]
end

function gradientOp(space::LagrangePolynomialSpace{<:Number,2})
    (nr, ns) = space.npoints
    (Dr, Ds) = space.deriv_mats

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()

    Dx = ⊗(Is, Dr)
    Dy = ⊗(Ds, Ir)

    DD = [Dx
          Dy]
end

function gradientOp(space::LagrangePolynomialSpace{<:Number,3})
    (Dr, Ds, Dt) = space.deriv_mats
    (nr, ns, nt) = space.npoints

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()
    It = IdentityOperator{nt}()

    Dx = ⊗(It, Is, Dr)
    Dy = ⊗(It, Ds, It)
    Dz = ⊗(Dt, Is, It)

    DD = [Dx
          Dy
          Dz]
end

###
# interpolation operators
###

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

    ⊗(Js, Jr)
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

    ⊗(Jt, Js, Jr)
end
#
