"""
Function Space Interface
"""
module Spaces

# TODO
#   - rename "Space" to Calculus
#   - figure out interplay between Space and Discretization
#       1. space is how you represent functions (with basis)
#       2. discretization is how you compute vector calculus operations

using SciMLOperators
import SciMLOperators: ⊗, IdentityOperator
import NNlib: gather, scatter
import FastGaussQuadrature: gausslobatto, gausslegendre, gausschebyshev
import FFTW: plan_rfft, plan_irfft
import UnPack: @unpack
import Plots

# overload
#import PDEInterfaces: dims, deform
import Base: eltype, length, summary, size
import Plots: plot, plot!

""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end

""" D-Dimensional spectral space """
abstract type AbstractSpectralSpace{T,D} <: AbstractSpace{T,D} end

""" D-Dimensional tensor-product space """
abstract type AbstractTensorProductSpace{T,D} <: AbstractSpace{T,D} end

include("interface.jl")
include("VectorCalculus.jl")
#include("GatherScatter.jl")

# Concrete Spaces
include("LagrangePolynomials/LagrangePolynomialSpace.jl")
include("TrigonometricPolynomials/Fourier.jl")

include("tensor.jl")
include("deform.jl")

export
       grid, domain,

       LagrangePolynomialSpace,
       GaussLobattoLegendre1D, GaussLegendre1D, GaussChebychev1D,
       GaussLobattoLegendre2D, GaussLegendre2D, GaussChebychev2D,

       gradOp, massOp, laplaceOp, advectionOp, divergenceOp,

       interpOp,

       deform

end
#
