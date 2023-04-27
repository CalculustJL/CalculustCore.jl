#
###
# Gather-Scatter Operators
###

"""
Q*Q'*u where Q: local -> global operator

TODO write GatherScatterOp that calls NNlib.gather, scatter
"""
function DSS(u, glo_num)
    Qu = NNlib.scatter(+, u, l2g) # Q
    QQtu = NNlib.gather(Qu, g2l)   # Q'

    return v
end

function gatherScatterOp(space::AbstractSpace)
    N = length(space)
    glo_num = global_numbering(space)

    # DSS
    op = (du, u, p, t) -> ()
    opi = (du, u, p, t) -> ()

    FunctionOperator(op; isinplace = true,
                     T = Bool,
                     size = (N, N), op_inverse = opi, issymmetric = true,
                     ishermitian = true)
end

function Qmatrix(n::Integer, periodic::Bool)
    Q = sparse(I, n, n - 1)

    if periodic
        Q[end, 1] = 1
    end

    Q
end

# replace with call to NNlib gather-scatter

#=
function GatherScatter(space::AbstractSpace)
    D = dims(space)
    N = length(space)

    domain = get_domain(space)
    periodic = isperiodic(domain)
    npoints = get_numpoints(space)

    if !prod(periodic...)
        return IdentityOp(N)
    end

    Qmats = Qmatrix.(npoints, periodic)

    Q = if D == dims(space)
        MatrixOp(Qmats...)
    else
        TensorProductOperator(Qmats...)

    QQt = Q * Q'
end
=#

#
