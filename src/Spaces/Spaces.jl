"""
Function Space Interface
"""
module Spaces

# TODO
#   - figure out interplay between Space and Discretization
#       1. space is how you represent functions (with basis)
#       2. discretization is how you compute vector calculus operations
#       3. have (space + discretization) --> calculus

using LinearAlgebra
using SciMLOperators
using SciMLOperators: IdentityOperator
import SparseArrays: sparse
using NNlib: gather, gather!, scatter, scatter!
using UnPack: @unpack
import Plots

using ..Domains
using ..Domains: AbstractDomain

# overload
import Base: eltype, length, size
import Base: summary, display, show
import Plots: plot, plot!
import ..Domains: dims, deform, ⊗

""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end

""" Spatial Discretizations """
abstract type AbstractDiscretization end

include("utils.jl")
include("interface.jl")
include("vectorcalculus.jl")
include("discretizations.jl")
include("gatherscatter.jl")

# Concrete Spaces
include("LagrangePolynomials/LagrangePolynomialSpace.jl")
include("TrigonometricPolynomials/Fourier.jl")

#include("tensor.jl") # TODO
include("transform.jl") # TODO
include("deform.jl")

export
       # interface
       dims, points, domain,
       local_numbering, global_numbering, basis, boundary_nodes,

       # discretizations
       Collocation, Galerkin,

       # operators
       massOp, gradOp, hessianOp, laplaceOp, advectionOp, divergenceOp,

       # interpolation
       interpOp,

       # lagrange polynomial space
       LagrangePolynomialSpace,
       GaussLobattoLegendre, GaussLegendre, GaussChebychev,

       # trigonometric polynomials
       FourierSpace,

       # misc
       deform,
       ⊗

end
#
