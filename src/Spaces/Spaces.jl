"""
Function Space Interface
"""
module Spaces

using Reexport

using LinearAlgebra
@reexport using SciMLOperators
using SciMLOperators: DEFAULT_UPDATE_FUNC
import SciMLOperators: IdentityOperator, NullOperator, ⊗

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
import ..Domains: dims, deform

""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end

""" Spatial Discretizations """
abstract type AbstractDiscretization end

include("utils.jl")
include("interface.jl")
include("vectorcalculus.jl")
include("discretizations.jl")
include("gatherscatter.jl")

#include("tensor.jl")
include("transform.jl")
include("deform.jl")

# Concrete Spaces
include("LagrangePolynomials/LagrangePolynomialSpace.jl")
include("TrigonometricPolynomials/Fourier.jl")

export
       ## Interface

       # from ..Domains
       dims, domain, deform,

       dims, domain, deform,
       # interface
       points, modes, basis, transform,
       local_numbering, global_numbering, boundary_nodes,
       # from SciMLOperators
       IdentityOperator, NullOperator, ⊗, cache_operator,

       ## Discretizations
       Collocation, Galerkin,

       ### Operators

       # vector calculus
       massOp, gradOp, hessianOp, laplaceOp, diffusionOp, advectionOp, divergenceOp,
       forcing,
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
