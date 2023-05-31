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
      bounding_box,
      isperiodic,
      boundaries,
      domain_tag,
      boundary_tags,
      ×,
      ⊗,
      deform,

      # concrete types
      IntervalDomain,

      UnitBoxDomain,
      UnitIntervalDomain,
      UnitSquareDomain,
      UnitCubeDomain,
      BoxDomain,

      ChebyshevDomain,
      FourierDomain,

      AnnulusDomain

end
#
