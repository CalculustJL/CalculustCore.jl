"""
Boundary Condition Interface
"""
module BoundaryConditions

#using Reexport
#@reexport using SciMLBase
import SciMLBase: init, solve; export init, solve

using SciMLOperators
using SciMLOperators: AbstractSciMLOperator, IdentityOperator
import UnPack: @unpack

using LinearAlgebra
using LinearSolve
using Plots

using ..Domains
using ..Domains: AbstractDomain

using ..Spaces
using ..Spaces: AbstractSpace

# overload
import SciMLBase: solve, init;
import Plots: plot, plot!

abstract type AbstractBoundaryCondition{T} end

abstract type AbstractBVPDEProblem <: SciMLBase.DEProblem end
abstract type AbstractBVPDECache <: SciMLBase.DECache end
abstract type AbstractBVPDEAlgorithm <: SciMLBase.DEAlgorithm end

include("utils.jl")

include("types.jl")
include("problem.jl")

export 
       # boundary conditions
       DirichletBC, NeumannBC, RobinBC, PeriodicBC,

       # boundary vale problem
       BVPDEProblem,

       # boundary value algorithms
       LinearBVPDEAlg, NonlinearBVPDEAlg

end
#
