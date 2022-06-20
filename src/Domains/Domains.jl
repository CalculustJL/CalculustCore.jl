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
#include("tensor.jl") # TODO - like annulus ⊗ interval
#include("meshed.jl") # TODO

include("deform.jl")
include("maps.jl")

export
       # interface
       dims,
       isperiodic,
       endpoints,
       boundary_tags,
       boundary_tag,
       num_boundaries,

       ⊗,
       ∈,
       in,

       # concrete types
       IntervalDomain,
       BoxDomain,

       PeriodicInterval,
       PeriodicBox,

       # conveniences
       reference_box,
       annulus_2D,

       deform

end
#
