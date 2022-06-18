#
"""
Deform domain, compute Jacobian of transformation, and its inverse

given 2D mapping

x = x(r,s), y = y(r,s)

compute

dXdR = [xr xs], dRdX = [rx ry], J = det(dXdR), Jinv = det(dRdX)
       [yr ys]         [sx sy]

"""
struct DeformedSpace{T,D,
                     Tspace<:AbstractSpace{T,D},
                     Tgrid,Tjacmat,Tjacimat,Tjac,Tjaci} <: AbstractSpace{T,D}
    space::Tspace
    grid::Tgrid # (x1, ..., xD,)
    dXdR::Tjacmat
    dRdX::Tjacimat
    J::Tjac
    Ji::Tjaci
    # add DeformedDomain
end

function deform(space::AbstractSpace{<:Number,D},
                mapping = nothing, isseparable = false
               ) where{D}
    if mappping === nothing
        J    = IdentityOp{D}()
        Jmat = Diagonal([J for i=1:D])
        @warn "mapping === nothing"
#       return space
        return DeformedSpace(space, grid(space), Jmat, Jmat, J, J)

    elseif isseparable # x = x(r), y = y(s)
        # eliminate cross terms by making dXdR, etc diagonal
    elseif rescaling # simple rescaling
        # make jac, dXdR, etc scaling operations
    end

    R = grid(space)
    X = mapping(R...)

    gradR = gradOp(space)

    """
    dXdR = [dx1/dr1 ... dx1/drD]
           [...     ...     ...]
           [dxD/dr1 ... dxD/drD]
    """
    dXdR = begin
        dXdR = gradR.(X) |> hcat
        dXdR = DiagonalOp.(dXdR)
        dXdR = dXdR'
    end

    J = if D == 1
        dXdR[1]
    elseif D == 2
        xr = dXdR[1]
        yr = dXdR[2]
        xs = dXdR[3]
        ys = dXdR[4]

        xr * ys - xs * yr
    elseif D == 3
        xr = dXdR[1]; yr = dXdR[2]; zr = dXdR[3]
        xs = dXdR[4]; ys = dXdR[5]; zs = dXdR[6]
        xt = dXdR[7]; yt = dXdR[8]; zt = dXdR[9]

        J = xr * (ys * zt - zs * yt) -
            xs * (yr * zt - zr * yt) +
            xt * (yr * zs - zr * ys)
    else
        det(dXdR) # errors
    end

    Ji = inv(D)

    dRdX = if D == 1
        fill(Ji, 1, 1)
    elseif D == 2
        xr = dXdR[1]
        yr = dXdR[2]
        xs = dXdR[3]
        ys = dXdR[4]

        rx =  (Ji * ys)
        ry = -(Ji * xs)
        sx = -(Ji * yr)
        sy =  (Ji * xr)

        dRdX = [rx ry 
                sx sy]
    elseif D == 3 # cramer's rule
        inv(dXdR)
    else
        inv(dXdR) # need even more robust AbstractOp interface
    end

    DeformedSpace(space, X, dXdR, dRdX, J, Ji)
end

get_grid(space::DeformedSpace) = space.grid

"""
[Dx] * u = [rx sx] * [Dr] * u
[Dy]     = [ry sy]   [Ds]
"""
function gradOp(space::DeformedSpace)
    gradR = gradOp(space.space)
    dRdX  = space.dRdX
    gradX = dRdX * gradR

    gradX
end

function massOp(space::DeformedSpace)
    M = massOp(space.space)
    J = space.J

    M * J
end

"""
(v,-∇² u) = (vx,ux) + (vy,uy)\n
         := a(v,u)\n
          = v' * A * u\n
          = (Q*R'*v)'*A_l*(Q*R'*u)\n
          = v'*R*Q'*A_l*Q*R'*u\n

implemented as

R'R * QQ' * A_l * u_loc
where A_l is

[Dr]'*[rx sx]'*[B 0]*[rx sx]*[Dr]
[Ds]  [ry sy]  [0 B] [ry sy] [Ds]

= [Dr]' * [G11 G12]' * [Dr]
  [Ds]    [G12 G22]    [Ds]
"""
function laplaceOp(space::DeformedSpace{<:Number, D}) where{D}

    Dr   = gradOp(space.space)
    M    = massOp(space)
    dRdX = space.dRdX

    MM = Diagonal([M for i=1:D])
    GG = if dRdX isa Diagonal
        dRdX' * MM * dRdX
    else
        dRdX' * MM * dRdX |> Symmetric # TODO avoid bottom 1/2 computation
    end

    laplOp = Dr' * GG * Dr

    return first(laplOp)
end

###
# Dealiased operators
###

function laplaceOp(space1::DeformedSpace{<:Number,D},
                   space2::DeformedSpace{<:Number,D},
                   J = nothing
                  ) where{D}
    J12 = J !== nothing ? J : interpOp(space1, space2)

    Dr1   = gradOp(space1.space)
    M2    = MassOp(space2)
    dRdX2 = space2.dRdX

    JD = J12 * Dr1

    MM2 = Diagonal([M2 for i=1:D])
    GG2 = dRdX' * MM2 * dRdX |> Symmetric

    laplOp = JD' * GG2 * JD

    return first(laplOp)
end
#
