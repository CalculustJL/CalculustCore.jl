#

struct BVPDEProblem{
                    isinplace,
                    T,
                    F,
                    fType,
                    uType,
                    Tbcs,
                    Tspace<:AbstractSpace{T},
                    P,
                    K,
                   } <: AbstractBVPDEProblem
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
    """Parameters"""
    p::P
    """Keyword arguments"""
    kwargs::K

    SciMLBase.@add_kwonly function BVPDEProblem(
         op::AbstractSciMLOperator,
         f,
         bc_dict::Dict,
         space::AbstractSpace{T},
         p = SciMLBase.NullParameters();
         u0 = nothing,
         kwargs...
        ) where{T}

        new{true,
            eltype(space),
            typeof(op),
            typeof(f),
            typeof(u0),
            typeof(bc_dict),
            typeof(space),
            typeof(p),
            typeof(kwargs)
           }(
             op, f, u0, bc_dict, space, p, kwargs
            )
    end
end

function Base.summary(io::IO, prob::BVPDEProblem)
    type_color, no_color = SciMLBase.get_colorizers(io)
    print(io,
          type_color, nameof(typeof(prob)),
          no_color," with uType ",
          type_color,typeof(prob.u0),
          no_color," with AType ",
          type_color,typeof(prob.op),
          no_color,". In-place: ",
          type_color,SciMLBase.isinplace(prob),
          no_color
         )
end

function Base.show(io::IO, mime::MIME"text/plain", A::BVPDEProblem)
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

struct BVPDECache{Top,Tu,Tbc,Tsp,Talg} <: AbstractBVPDECache
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
    """ Algorithm """
    alg::Talg
end

struct BVPDESolution{T,D,uType,R,A,C} #<: SciMLBase.AbstractDAESolution{T,D}
    u::uType
    resid::R
    alg::A
    retcode::Symbol
    iters::Int
    cache::C

    function BVPDESolution(u, resid, alg, retcode, iters, cache)
        @unpack space = cache

        new{
            eltype(space),
            dims(space),
            typeof(u),
            typeof(resid),
            typeof(alg),
            typeof(cache),
           }(
             u, resid, alg, retcode, iters, cache,
            )
    end
end

function build_bvpde_solution(alg, u, resid, cache; retcode=:Default, iters=0)
    BVPDESolution(u, resid, alg, retcode, iters, cache)
end

function Plots.plot(sol::BVPDESolution{<:Number,1})
    plot(sol.u, sol.cache.space)
end

function Plots.plot(sol::BVPDESolution{<:Number,2}; a=45, b=60)
    plot(sol.u, sol.cache.space)
end

Base.@kwdef struct LinearBVPDEAlg{Tl} <: AbstractBVPDEAlgorithm
    linalg::Tl = nothing
end

#TODO integrate NonlinearSolve.jl with LinearSolve.jl first
Base.@kwdef struct NonlinearBVPDEAlg{Tnl} <: AbstractBVPDEAlgorithm
    nlalg::Tnl = nothing
end

function makeLHS(op::AbstractSciMLOperator, bc::AbstractBoundaryCondition)
    @unpack mask_dir, amask_dir = bc

    #TODO
    """
    how to empty dirichlet row? -
        https://github.com/vpuri3/PDEInterfaces.jl/issues/1

    is having the construct u := u_inhom + u_hom the only way?
    then would have to interpolate u_inhom into interior.
    what are the consequences?
    """

    lhs = mask_dir * op + amask_dir
end

function makeRHS(f, bc::AbstractBoundaryCondition)
    @unpack bc_dict, antimasks, mask_dir, space = bc

    M = massOp(space)
    b = M * f

    dirichlet = zero(b)
    neumann   = zero(b)
    robin     = zero(b)

    pts = grid(space)
    dom = domain(space)

    for i=1:num_boundaries(dom)
        tag   = boundary_tag(dom, i)
        bc    = bc_dict[tag]
        amask = antimasks[i]

        if bc isa DirichletBC
            dirichlet += amask * bc.f(pts...)
        elseif bc isa NeumannBC
            neumann += amask * bc.f(pts...)
        elseif bc isa RobinBC
            # TODO
        elseif bc isa PeriodicBC
            continue
        end
    end

    b = (mask_dir * b) + dirichlet - neumann + robin
end

function SciMLBase.solve(cache::BVPDECache; kwargs...)
    @unpack op, f, u, bc, space, alg = cache
    @unpack linalg = alg

    lhsOp = makeLHS(op, bc)
    lhsOp = cache_operator(lhsOp, f)
    rhs   = makeRHS(f, bc)

    linprob = LinearProblem(lhsOp, rhs; u0=_vec(u))
    linsol  = solve(linprob, linalg; kwargs...)

    resid = norm(lhsOp * linsol.u - rhs, Inf)

    build_bvpde_solution(alg, u, resid, cache; iters=linsol.iters)
end

function SciMLBase.init(prob::AbstractBVPDEProblem,
                        alg::AbstractBVPDEAlgorithm = nothing;
                        abstol=default_tol(eltype(prob.op)),
                        reltol=default_tol(eltype(prob.op)),
                        maxiters=length(prob.f),
                        verbose=false,
                        kwargs...
                       )
    @unpack op, f, u0, bc_dict, space = prob

    alg = alg isa Nothing ? LinearBVPDEAlg() : alg

    u  = u0 isa Nothing ? zero(f) : u0
    bc = BoundaryCondition(bc_dict, space)

    BVPDECache(op, f, u, bc, space, alg)
end

function SciMLBase.solve(prob::BVPDEProblem, args...; kwargs...)
    solve(init(prob, nothing, args...; kwargs...))
end

function SciMLBase.solve(prob::BVPDEProblem, alg::Union{AbstractBVPDEAlgorithm,Nothing}, args...; kwargs...)
    solve(init(prob, alg, args...; kwargs...); kwargs...)
end
#
