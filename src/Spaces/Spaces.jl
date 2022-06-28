"""
Function Space Interface
"""
module Spaces

using Reexport
using LinearAlgebra
using Plots: Plots

@reexport using SciMLOperators
using SciMLOperators: AbstractSciMLOperator, DEFAULT_UPDATE_FUNC

using UnPack: @unpack
using Setfield: @set!
using NNlib: gather, gather!, scatter, scatter!
import SparseArrays: sparse
import FFTW
import CUDA

using ..Domains
using ..Domains: AbstractDomain

# overload
import Base: eltype, length, size
import Base: summary, display, show
import Plots: plot, plot!
import ..Domains: dims, deform
import SciMLOperators: IdentityOperator, NullOperator, ⊗

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
       ### Interface
       points,
       modes,
       basis,
       local_numbering,
       global_numbering,
       boundary_nodes,

       transform,
       make_transform,

       # from ..Domains
       dims,
       domain,
       deform,

       # from SciMLOperators
       IdentityOperator,
       NullOperator,
       ⊗,

       ### Discretizations
       Collocation,
       Galerkin,

       ### Operators

       # vector calculus
       massOp,
       gradientOp,
       hessianOp,
       laplaceOp,
       biharmonicOp,
       diffusionOp,
       advectionOp,
       divergenceOp,
       forcingOp,

       # interpolation
       interpOp,

       # transform
       transformOp,
       truncationOp,

       ### Spaces

       # Lagrange polynomial spaces
       LagrangePolynomialSpace,
       GaussLobattoLegendre,
       GaussLegendre,
       GaussChebychev,

       # Trigonometric polynomial spaces
#      SineSpace,
#      CosineSpace,
       FourierSpace

end
#
