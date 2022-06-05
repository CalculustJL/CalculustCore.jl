#
###
# Gather-Scatter Operators - enforce continuity/ periodicity
###

abstract type AbstractGatherScatterOperator{D} <: SciMLOperators.AbstractSciMLOperator{Bool} end

Base.adjoint(A::AbstractGatherScatterOperator) = A

# TODO write GatherScatterOp that calls NNlib.gather, scatter
# implement mul!, *
struct GatherScatter{D} <: AbstractGatherScatterOperator{D}
    global_numbering
    implementation
end

"""
Q*Q'*u where Q: local -> global operator
"""
function DSS(u,l2g,g2l)

    Qu   = NNlib.scatter(+,u,l2g) # Q
    QQtu = NNlib.gather(Qu,g2l)   # Q'

    return v
end

"""
map from local vector to global vector
"""
function Qmatrix(n::Integer, periodic::Bool)
    Q = sparse(I,n, n-1)

    if periodic
        Q[end,1] = 1
    end

    Q
end

function GatherScatter(space::AbstractSpectralSpace{<:Number,D}) where{D}
    domain = get_domain(space)
    periodic = isperiodic(domain)
    npoints = get_numpoints(space)

    if !prod(periodic...)
        return IdentityOp{D}()
    end

    Qmats = Qmatrix.(npoints, periodic)

    Q = if D == 1
        MatrixOp(Qmats...)
    elseif D == 2
        TensorProductOp2D(Qmats...)
    elseif D == 3
        TensorProductOp3D(Qmats...)
    end

    QQt = Q * Q' # replace with call to NNlib gather-scatter
end
#
