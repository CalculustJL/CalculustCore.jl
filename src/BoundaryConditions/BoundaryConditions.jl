"""
Boundary Condition Interface
"""
module BoundaryConditions # change to BoundaryValueProblems

using DocStringExtensions

using Reexport
@reexport using SciMLBase
using SciMLOperators: AbstractSciMLOperator, IdentityOperator

import UnPack: @unpack

using LinearAlgebra
using LinearSolve

using ..Domains
using ..Domains: AbstractDomain

using ..Spaces
using ..Spaces: AbstractSpace, AbstractDiscretization

"""
$TYPEDEF
"""
abstract type AbstractBoundaryCondition{T} end

"""
$TYPEDEF
"""
abstract type AbstractBoundaryValueProblem <: SciMLBase.DEProblem end

"""
$TYPEDEF
"""
abstract type AbstractBoundaryValueCache <: SciMLBase.DECache end

"""
$TYPEDEF
"""
abstract type AbstractBoundaryValueAlgorithm <: SciMLBase.DEAlgorithm end

include("utils.jl")
include("conditions.jl")

include("types.jl")
include("make_lhs_rhs.jl")
include("problem.jl")

export
      # boundary conditions
      DirichletBC, NeumannBC, RobinBC, PeriodicBC,

      # boundary vale problem
      BoundaryValueProblem,

      # boundary value algorithms
      LinearBoundaryValueAlg, NonlinearBoundaryValueAlg

end
#
