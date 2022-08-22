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
using Lux: cpu, gpu
import SparseArrays: sparse

using ..Domains
using ..Domains: AbstractDomain

# overload
import Base: eltype, length, size
import Base: summary, display, show

import SciMLOperators: IdentityOperator, NullOperator, ⊗
import ..Domains: dims, deform

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

include("NDgrid.jl") # TODO - use https://github.com/JuliaArrays/LazyGrids.jl instead

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

       ### from Lux
       cpu,
       gpu,

       ### Discretizations
       Collocation,
       Galerkin,

       ### Operators

       # from SciMLOperators
       IdentityOperator,
       NullOperator,
       ⊗,

       ### operators

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

       transformOp,
       truncationOp

end
#
