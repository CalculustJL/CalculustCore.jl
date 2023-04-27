"""
Boundary Condition Interface
"""
module BoundaryConditions # change to BoundaryValueProblems

#using Reexport
#@reexport using SciMLBase
import SciMLBase: init, solve
export init, solve

using SciMLOperators
using SciMLOperators: AbstractSciMLOperator, IdentityOperator
import UnPack: @unpack

using LinearAlgebra
using LinearSolve
using Plots

using ..Domains
using ..Domains: AbstractDomain

using ..Spaces
using ..Spaces: AbstractSpace, AbstractDiscretization

# overload
import SciMLBase: solve, init
import Plots: plot, plot!

abstract type AbstractBoundaryCondition{T} end

abstract type AbstractBoundaryValueProblem <: SciMLBase.DEProblem end
abstract type AbstractBoundaryValueCache <: SciMLBase.DECache end
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
