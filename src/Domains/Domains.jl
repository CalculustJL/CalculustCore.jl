"""
Abstract Domain Interface
"""
module Domains

using DocStringExtensions
using Setfield: @set!

import Base: show
import Base: ndims, eltype, ∈, in
import LinearAlgebra: ×

"""
$TYPEDEF

`D`-dimensional domain of type `T`
"""
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
      expanse,
      bounding_box,
      isperiodic,
      boundaries,
      domain_tag,
      boundary_tags,
      num_boundaries,
      ×,
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
