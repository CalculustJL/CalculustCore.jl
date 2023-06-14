"""
Abstract interface for function spaces
"""
module Spaces

using DocStringExtensions

using Reexport
using LinearAlgebra

# operator algebra
@reexport using SciMLOperators
using SciMLOperators: AbstractSciMLOperator, DEFAULT_UPDATE_FUNC,
                      IdentityOperator, NullOperator, ⊗

# caching
using Setfield: @set!

# gather-scatter
# using NNlib: gather, gather!, scatter, scatter!
# import SparseArrays: sparse

@reexport using ..Domains

# interface
import Base: show
import Base: ndims, eltype, length, size

"""
$TYPEDEF

Function space over a `D`-Dimensional domain
"""
abstract type AbstractSpace{T, D} end

"""
$TYPEDEF

Spatial discretizations scheme
"""
abstract type AbstractDiscretization end

include("utils.jl")

# AbstractSpace interface
include("interface.jl")

# PDE operators
include("vectorcalculus.jl")

include("discretizations/collocation.jl")
include("discretizations/galerkin.jl")

include("NDgrid.jl")
# include("gatherscatter.jl")

# include("tensor.jl")
include("transform.jl")
include("deform.jl")

export
      ### from ..Domains
      deform,

      ### Interface
      domain,
      points,
      global_numbering,
      boundary_nodes, ndgrid,

      modes,
      mode_size,
      transform,
      make_transform,

      # SciMLOperators
      ⊗,
      IdentityOperator,
      NullOperator,

      ### Discretizations
      Collocation,
      Galerkin,

      ### Operators
      IdentityOperator,
      NullOperator, massOp,
      gradientOp,
      hessianOp,
      laplaceOp,
      biharmonicOp,
      diffusionOp,
      advectionOp,
      divergenceOp,
      forcingOp, interpOp, transformOp,
      truncationOp

end
#
