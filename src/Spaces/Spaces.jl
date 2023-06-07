"""
Abstract interface for function spaces
"""
module Spaces

using Reexport
using LinearAlgebra

# operator algebra
@reexport using SciMLOperators
using SciMLOperators: AbstractSciMLOperator, DEFAULT_UPDATE_FUNC,
                      IdentityOperator, NullOperator

# caching
using UnPack: @unpack
using Setfield: @set!

# gather-scatter
# using NNlib: gather, gather!, scatter, scatter!
import SparseArrays: sparse

@reexport using ..Domains

# interface
import Base: summary, show
import Base: eltype, length, size

"""
Function space over a `D`-Dimensional domain
"""
abstract type AbstractSpace{T, D} end

"""
Spatial Discretizations
"""
abstract type AbstractDiscretization end

include("utils.jl")

# interface
include("interface.jl")

# operators
include("vectorcalculus.jl")
include("discretizations.jl")
include("gatherscatter.jl")

include("NDgrid.jl") # TODO - use LazyGrids.jl instead

#include("tensor.jl")
include("transform.jl")
include("deform.jl")

export
      ### from ..Domains
      dims,
      deform,

      ### Interface
      points,
      modes,
      mode_size,
      #basis,
      domain,
      mass_matrix,
      global_numbering,
      boundary_nodes, ndgrid, transform,
      make_transform,

      # from SciMLOperators
      ⊗,

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
