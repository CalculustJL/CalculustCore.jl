"""
Function Space Interface
"""
module Spaces

using Reexport

using LinearAlgebra
@reexport using SciMLOperators
using SciMLOperators: IdentityOperator, NullOperator
import SparseArrays: sparse
using NNlib: gather, gather!, scatter, scatter!
using UnPack: @unpack
using Lazy: @forward
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
       ## Interface

       # from ..Domains
       dims, domain, deform, ⊗,
       # interface
       points, modes, basis, transform,
       local_numbering, global_numbering, boundary_nodes,
       # from SciMLOperators
       cache_operator,

       ## Discretizations
       Collocation, Galerkin,

       ### Operators

       # vector calculus
       massOp, gradOp, hessianOp, laplaceOp, diffusionOp, advectionOp, divergenceOp,
       # interpolation
       interpOp,

       ### Spaces

       # Lagrange polynomial spaces
       LagrangePolynomialSpace,
       GaussLobattoLegendre, GaussLegendre, GaussChebychev,
       # Trigonometric polynomial spaces
       FourierSpace #, SineSpace, CosineSpace,

end
#
