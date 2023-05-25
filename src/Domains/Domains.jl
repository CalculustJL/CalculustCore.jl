"""
Domain Interface
"""
module Domains

using DocStringExtensions
using Setfield

import Base: eltype, ∈, in
import SciMLOperators: ⊗
import LinearAlgebra: ×

""" D-Dimensional Domain """
abstract type AbstractDomain{T<:Number, D} end

# Interface
include("interface.jl")

# Concrete Types
include("basic.jl")
include("product.jl")
#include("meshed.jl") # TODO

include("deform.jl")
include("concrete.jl")

export

      # interface
      dims,
      expanse,
      isperiodic,
      boundaries,
      domain_tag,
      boundary_tag,
      num_boundaries,
      bounding_box,
      ×, ⊗,
      deform,

      # concrete types
      PointDomain,
      IntervalDomain,
      ChebyshevDomain,
      FourierDomain,
      AnnulusDomain

end
#
