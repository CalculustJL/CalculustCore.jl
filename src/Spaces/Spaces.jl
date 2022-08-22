"""
Function Space Interface
"""
module Spaces

using Reexport
using LinearAlgebra
using Plots: Plots

# operator algebra
@reexport using SciMLOperators
using SciMLOperators: AbstractSciMLOperator, DEFAULT_UPDATE_FUNC

# caching
using UnPack: @unpack
using Setfield: @set!

# GPU
using Lux: cpu, gpu

# gather-scatter
using NNlib: gather, gather!, scatter, scatter!
import SparseArrays: sparse

using ..Domains

# interface
import Base: eltype, length, size
import Base: summary, display, show

import ..Domains: dims, deform
import SciMLOperators: IdentityOperator, NullOperator, ⊗

# plot recipies
import Plots: plot, plot!, @animate, animate

""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end

""" Spatial Discretizations """
abstract type AbstractDiscretization end

include("utils.jl")
include("interface.jl")
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
       domain,
       deform,

       ### Interface
       points,
       modes,
       basis,
       local_numbering,
       global_numbering,
       boundary_nodes,

       ndgrid,

       transform,
       make_transform,

       # from SciMLOperators
       ⊗,

       ### from Lux
       cpu,
       gpu,

       ### Discretizations
       Collocation,
       Galerkin,

       ### Operators
       IdentityOperator,
       NullOperator,

       massOp,
       gradientOp,
       hessianOp,
       laplaceOp,
       biharmonicOp,
       diffusionOp,
       advectionOp,
       divergenceOp,
       forcingOp,

       interpOp,

       transformOp,
       truncationOp

end
#
