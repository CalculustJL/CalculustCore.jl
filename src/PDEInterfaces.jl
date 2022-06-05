module PDEInterfaces

using Reexport

##########################################
#@reexport using SciMLBase
import SciMLBase
import SciMLBase:solve, init
export solve, init

using SciMLOperators
import SciMLOperators: ⊗, IdentityOperator
import SciMLOperators: AbstractSciMLOperator, _reshape, _vec
##########################################

using LinearAlgebra
using LinearSolve

import UnPack: @unpack
import Setfield: @set!
import SparseArrays: sparse
import NNlib: gather, scatter
import FastGaussQuadrature: gausslobatto, gausslegendre, gausschebyshev
import FFTW: plan_rfft, plan_irfft

###
# Abstract Supertypes
###

""" D-Dimensional physical domain """
abstract type AbstractDomain{T,D} end
""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end

###
# AbstractSpace subtypes
###

""" D-Dimensional spectral space """
abstract type AbstractSpectralSpace{T,D} <: AbstractSpace{T,D} end

""" D-Dimensional tensor-product space """
abstract type AbstractTensorProductSpace{T,D} <: AbstractSpace{T,D} end

AbstractSupertypes{T,D} = Union{
                                AbstractSpace{T,D},
                                AbstractDomain{T,D},
                               }

# traits
dims(::AbstractSupertypes{T,D}) where{T,D} = D
Base.eltype(::Union{
                    AbstractSpace{T,D},
                    AbstractDomain{T,D},
                   }
           ) where{T,D} = T

# misc
include("NDgrid.jl")

# domain
include("Domain.jl")
include("DomainMaps.jl")

# misc
include("BoundaryConditions.jl")
#include("GatherScatter.jl")

# space
include("Space.jl")
include("DeformSpace.jl")

# discretizations
#include("Discretizations.jl")

# polynomial spaces
include("LagrangeMatrices.jl")
include("LagrangePolynomialSpace.jl")

# fourier spaces
#include("FourierSpace.jl")

# problems
include("BoundaryValueProblem.jl")
#include("EigenValueProblem.jl")

export 
       dims,

       ## Domains
       isperiodic, endpoints, boundary_tags, boundary_tag,
       num_boundaries, bounding_box,

       IntervalDomain, BoxDomain, deform,
       reference_box, annulus_2D,

       ## Fields
       Field,

       ## Spaces
       get_grid, get_domain, numpoints,

       gradOp, massOp, laplaceOp, advectionOp, divergenceOp,

       interpOp,

       LagrangePolynomialSpace,
       GaussLobattoLegendre1D, GaussLegendre1D, GaussChebychev1D,
       GaussLobattoLegendre2D, GaussLegendre2D, GaussChebychev2D,

       ## Boundary conditions
       DirichletBC, NeumannBC, RobinBC,

       ## Boundary vale problem
       BoundaryValuePDEProblem, LinearBVPDEAlg, NonlinearBVPDEAlg

end # module
