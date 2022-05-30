module PDEInterfaces

using Reexport

@reexport using SciMLBase

using LinearAlgebra
using LinearSolve

import UnPack: @unpack
import Setfield: @set!
import Base.ReshapedArray
import SciMLBase: AbstractDiffEqOperator
import Lazy: @forward
import SparseArrays: sparse
import NNlib: gather, scatter
import FastGaussQuadrature: gausslobatto, gausslegendre, gausschebyshev
import FFTW: plan_rfft, plan_irfft

# AbstractVector subtyping
import Base: summary, show, similar, zero
import Base: size, getindex, setindex!, IndexStyle
import Base.Broadcast: BroadcastStyle

# overload maths
import Base: +, -, *, /, \, adjoint, ∘, inv, one, convert
import LinearAlgebra: mul!, ldiv!, lmul!, rmul!

###
# Abstract Supertypes
###

""" Scalar function field in D-Dimensional space """
abstract type AbstractField{T,D} <: AbstractVector{T} end
""" Operators acting on fields in D-Dimensional space """
abstract type AbstractOperator{T,D} <: AbstractDiffEqOperator{T} end
""" D-Dimensional physical domain """
abstract type AbstractDomain{T,D} end
""" Function space in D-Dimensional space """
abstract type AbstractSpace{T,D} end
""" Boundary condition on domain in D-Dimensional space """
abstract type AbstractBonudaryCondition{T,D} end

###
# AbstractOperator subtypes
###

""" Tensor product operator in D-Dimensional space """
abstract type AbstractTensorProductOperator{T,D} <: AbstractOperator{T,D} end

""" Gather-Scatter operator in D-Dimensional space """
abstract type AbstractGatherScatterOperator{T,D} <: AbstractOperator{T,D} end

###
# AbstractSpace subtypes
###

""" D-Dimensional spectral space """
abstract type AbstractSpectralSpace{T,D} <: AbstractSpace{T,D} end

""" D-Dimensional tensor-product space """
abstract type AbstractTensorProductSpace{T,D} <: AbstractSpace{T,D} end

###
# other abstract types
###

""" Boundary Condition on D-Dimensional domain """
abstract type AbstractBoundaryCondition{T,D} end

AbstractSupertypes{T,D} = Union{
                                AbstractField{T,D},
                                AbstractOperator{T,D},
                                AbstractSpace{T,D},
                                AbstractDomain{T,D},
                                AbstractBonudaryCondition{T,D}
                               }

# traits

dims(::AbstractSupertypes{T,D}) where{T,D} = D
Base.eltype(::Union{
                    AbstractSpace{T,D},
                    AbstractDomain{T,D},
                    AbstractBoundaryCondition{T,D},
                   }
           ) where{T,D} = T

# misc
include("utils.jl")
include("NDgrid.jl")

# field
include("Field.jl")

# operators
include("OperatorBasics.jl")
include("Operators.jl")

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
