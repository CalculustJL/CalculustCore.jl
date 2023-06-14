#
struct BoundaryValueProblem{
                            isinplace,
                            T,
                            F,
                            fType,
                            uType,
                            Tbcs,
                            Tspace <: AbstractSpace{T},
                            Tdiscr <: AbstractDiscretization,
                            P,
                            K
                            } <: AbstractBoundaryValueProblem
    """Neumann Operator"""
    op::F
    """Right-hand-side forcing vector"""
    f::fType
    """Initial guess"""
    u0::uType
    """Boundary condition object"""
    bc_dict::Tbcs
    """Function space"""
    space::Tspace
    """Discretization"""
    discr::Tdiscr
    """Parameters"""
    p::P
    """Keyword arguments"""
    kwargs::K

    SciMLBase.@add_kwonly function BoundaryValueProblem(op::AbstractSciMLOperator,
                                                        f,
                                                        bc_dict::Dict,
                                                        space::AbstractSpace{T},
                                                        discr::AbstractDiscretization,
                                                        p = SciMLBase.NullParameters();
                                                        u0 = nothing,
                                                        kwargs...) where {T}
        new{true,
            eltype(space),
            typeof(op),
            typeof(f),
            typeof(u0),
            typeof(bc_dict),
            typeof(space),
            typeof(discr),
            typeof(p),
            typeof(kwargs)
            }(op, f, u0, bc_dict, space, discr, p, kwargs)
    end
end

function Base.summary(io::IO, prob::BoundaryValueProblem)
    type_color, no_color = SciMLBase.get_colorizers(io)
    print(io,
          type_color, nameof(typeof(prob)),
          no_color, " with uType ",
          type_color, typeof(prob.u0),
          no_color, " with AType ",
          type_color, typeof(prob.op),
          no_color, ". In-place: ",
          type_color, SciMLBase.isinplace(prob),
          no_color)
end

function Base.show(io::IO, mime::MIME"text/plain", A::BoundaryValueProblem)
    summary(io, A)
    println(io)
    println(io, "u0: ")
    show(io, mime, A.u0)
    println(io, "Neumann Operator: ")
    show(io, mime, A.op)
    println(io, "Forcing Vector: ")
    show(io, mime, A.f)
    println(io, "Boundary Conditions: ")
    show(io, mime, A.bc_dict)
    println(io, "Function Space: ")
    show(io, mime, A.space)
end

struct BoundaryValueCache{Top, Tu, Tbc, Tsp, Tdi, Talg} <: AbstractBoundaryValueCache
    """Neumann Operator"""
    op::Top
    """Right-hand-side forcing vector"""
    f::Tu
    """Initial guess"""
    u::Tu
    """Boundary condition object"""
    bc::Tbc
    """Function space"""
    space::Tsp
    """Discretization"""
    discr::Tdi
    """Algorithm"""
    alg::Talg
end

struct BoundaryValueSolution{T, D, uType, R, A, C} #<: SciMLBase.AbstractDAESolution{T,D}
    u::uType
    resid::R
    alg::A
    retcode::Symbol
    iters::Int
    cache::C

    function BoundaryValueSolution(u, resid, alg, retcode, iters, cache)
        @unpack space = cache

        new{
            eltype(space),
            ndims(space),
            typeof(u),
            typeof(resid),
            typeof(alg),
            typeof(cache)
            }(u, resid, alg, retcode, iters, cache)
    end
end

function build_bv_solution(alg, u, resid, cache; retcode = :Default, iters = nothing)
    BoundaryValueSolution(u, resid, alg, retcode, iters, cache)
end

Base.@kwdef struct LinearBoundaryValueAlg{Tl} <: AbstractBoundaryValueAlgorithm
    linalg::Tl = nothing
end

#TODO integrate NonlinearSolve.jl with LinearSolve.jl first
Base.@kwdef struct NonlinearBoundaryValueAlg{Tnl} <: AbstractBoundaryValueAlgorithm
    nlalg::Tnl = nothing
end

function SciMLBase.solve(cache::BoundaryValueCache; kwargs...)
    @unpack op, f, u, bc, space, alg = cache
    @unpack linalg = alg

    lhsOp = makeLHS(op, bc)
    lhsOp = cache_operator(lhsOp, f)
    rhs = makeRHS(f, bc)

    linprob = LinearProblem(lhsOp, rhs; u0 = vec(u))
    linsol = solve(linprob, linalg; kwargs...)

    resid = norm(lhsOp * linsol.u - rhs, Inf)

    build_bv_solution(alg, u, resid, cache; iters = linsol.iters)
end

function SciMLBase.init(prob::AbstractBoundaryValueProblem,
                        alg::AbstractBoundaryValueAlgorithm = nothing;
                        abstol = default_tol(eltype(prob.op)),
                        reltol = default_tol(eltype(prob.op)),
                        maxiters = length(prob.f),
                        verbose = false,
                        kwargs...)
    @unpack op, f, u0, bc_dict, space, discr = prob

    alg = alg isa Nothing ? LinearBoundaryValueAlg() : alg

    u = u0 isa Nothing ? zero(f) : u0
    bc = BoundaryCondition(bc_dict, space, discr)

    BoundaryValueCache(op, f, u, bc, space, discr, alg)
end

function SciMLBase.solve(prob::BoundaryValueProblem, args...; kwargs...)
    solve(init(prob, nothing, args...; kwargs...))
end

function SciMLBase.solve(prob::BoundaryValueProblem,
                         alg::Union{AbstractBoundaryValueAlgorithm, Nothing}, args...;
                         kwargs...)
    solve(init(prob, alg, args...; kwargs...); kwargs...)
end
#
