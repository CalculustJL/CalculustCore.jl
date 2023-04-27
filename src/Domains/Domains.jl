"""
Domain Interface
"""
module Domains

import Base: eltype, ∈, in
import SciMLOperators: ⊗

""" D-Dimensional Domain """
abstract type AbstractDomain{T, D} end

# Interface
include("interface.jl")

#include("tensor.jl") # TODO - like annulus ⊗ interval. for extruding

# Concrete Types
include("interval.jl")
include("box.jl")
#include("meshed.jl") # TODO

include("deform.jl")
include("maps.jl")

export

      # interface
      dims,
      lengths,
      isperiodic,
      endpoints,
      boundary_tags,
      boundary_tag,
      num_boundaries,
      bounding_box, ⊗,
      deform,

      # concrete types
      IntervalDomain,
      BoxDomain,
      GaussLobattoLegendreDomain,
      ChebyshevDomain,
      FourierDomain,
      AnnulusDomain

end
#
