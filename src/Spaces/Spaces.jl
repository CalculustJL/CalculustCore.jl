"""
Function Space Interface
"""
module Spaces # replace with calculus (?)
#
# TODO
#   - rename "Space" to Calculus
#   - figure out interplay between Space and Discretization
#       1. space is how you represent functions (with basis)
#       2. discretization is how you compute vector calculus operations

using SciMLOperators
import SciMLOperators: ⊗, IdentityOperator

import PDEInterfaces: dims, deform
import Base: eltype, length, summary, size
import Plots: plot

""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end

# Interface
include("interface.jl")

# Concrete Spaces
include("LagrangePolynomials/LagrangePolynomialSpace.jl")
#include("Trig/Fourier.jl")

include("tensor.jl")
include("deform.jl")

end
#
