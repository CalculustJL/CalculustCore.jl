"""
Domain Interface
"""
module Domains

import Base: eltype, ∈, in
import SciMLOperators: ⊗

""" D-Dimensional Domain """
abstract type AbstractDomain{T,D} end

# Interface
include("interface.jl")

# Concrete Types
include("interval.jl")
include("box.jl")
#include("tensor.jl") # TODO - like annulus ⊗ interval. for extruding
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
       bounding_box,

       ⊗,
       ∈,
       in,

       # concrete types
       IntervalDomain,
       BoxDomain,

       # conveniences
       GaussLobattoLegendreDomain, ChebychevDomain, FourierDomain, AnnulusDomain,

       deform

end
#
