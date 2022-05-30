#
abstract type AbstractBoundaryValuePDEProblem <: SciMLBase.DEProblem end
abstract type AbstractBoundaryValuePDECache <: SciMLBase.DECache end
abstract type AbstractBoundaryValuePDEAlgorithm <: SciMLBase.DEAlgorithm end

struct BoundaryValuePDEProblem{
                               isinplace,F,fType,uType,Tbcs,Tspace,P,K,
                              } <: AbstractBoundaryValuePDEProblem
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

    SciMLBase.@add_kwonly function BoundaryValuePDEProblem(
         op::AbstractOperator{<:Number,D},
         f::AbstractField{<:Number,D},
         bc_dict::Dict,
         space::AbstractSpace{<:Number,D},
         p = SciMLBase.NullParameters();
         u0::Union{AbstractField{<:Number,D},Nothing} = nothing,
         kwargs...
        ) where{D}

        new{true,
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

function Base.summary(io::IO, prob::AbstractBoundaryValuePDEProblem)
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

function Base.show(io::IO, mime::MIME"text/plain", A::AbstractBoundaryValuePDEProblem)
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

struct BoundaryValuePDECache{Top,Tu,Tbc,Tsp,Talg} <: AbstractBoundaryValuePDECache
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

Base.@kwdef struct LinearBVPDEAlg{Tl} <: AbstractBoundaryValuePDEAlgorithm
    linalg::Tl = nothing
end

#TODO integrate NonlinearSolve.jl with LinearSolve.jl first
Base.@kwdef struct NonlinearBVPDEAlg{Tnl} <: AbstractBoundaryValuePDEAlgorithm
    nlalg::Tnl = nothing
end

function makeLHS(op::AbstractOperator{<:Number,D},
                 bc::AbstractBoundaryCondition{<:Number,D}) where{D}
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

function makeRHS(f::AbstractField{<:Number,D},
                 bc::AbstractBoundaryCondition{<:Number,D}) where{D}
    @unpack bc_dict, antimasks, dirichlet_mask, space = bc

    M = MassOp(space)
    b = M * f

    dirichlet = zero(b)
    neumann   = zero(b)
    robin     = zero(b)

    grid = get_grid(space)

    for i=1:num_boundaries(domain)
        tag   = boundary_tag(domain, i)
        bc    = bc_dict[tag]
        amask = antimasks[i]

        if bc isa DirichletBC
            dirichlet += amask * bc.f(grid...)
        elseif bc isa NeumannBC
            neumann += amask * bc.f(grid...)
        elseif bc isa RobinBC
            # TODO
        elseif bc isa PeriodicBC
            continue
        end
    end

    b = dirichlet_mask * b + dirichlet - neumann + robin
end

function SciMLBase.solve(cache::BoundaryValuePDECache)
    @unpack op, f, u, bc, space, alg = cache
    @unpack linalg = alg

    lhsOp = makeLHS(op, bc)
    rhs   = makeRHS(f, bc)

    linprob = LinearProblem(lhsOp, b; u0=u)
    linsol  = solve(linprob, linalg)

    linsol.u
end

function SciMLBase.init(prob::AbstractBoundaryValuePDEProblem, alg::AbstractBoundaryValuePDEAlgorithm = nothing)
    @unpack op, f, u0, bc_dict, space = prob

    u  = u0 isa Nothing ? zero(f) : u0
    bc = BoundaryCondition(bc_dict, space)

    alg = alg isa Nothing ? LinearBVPDEAlg() : alg

    BoundaryValuePDECache(op, f, u, bc, space, alg)
end

function SciMLBase.solve(prob::BoundaryValuePDEProblem, alg::AbstractBoundaryValuePDEAlgorithm)
    solve(init(prob, alg))
end

function SciMLBase.solve(cache::BoundaryValuePDECache, alg::AbstractBoundaryValuePDEAlgorithm)
    solve(cache, cache.alg)
end
#
