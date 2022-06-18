module PDEInterfaces

using Reexport

#@reexport using SciMLBase
import SciMLBase; import SciMLBase:solve, init; export solve, init

using LinearAlgebra
using LinearSolve
import Plots

import Plots: plot, plot!
import UnPack: @unpack

# traits
dims(::AbstractArray{<:Any,D}) where{D} = D

# misc
include("Utils/utils.jl")

include("Domains/Domains.jl")
include("Spaces/Spaces.jl")

#include("Field.jl")

include("BoundaryConditions.jl")
#include("Discretizations.jl")

# problems
include("BoundaryValueProblem.jl")
#include("EigenValueProblem.jl")

export 
       dims,

       ## Plots
       plot, plot!

       ## Fields
#      Field,

       ## Spaces
       grid, domain,

       ## Boundary conditions
       DirichletBC, NeumannBC, RobinBC, PeriodicBC,

       ## Boundary vale problem
       BVPDEProblem, LinearBVPDEAlg, NonlinearBVPDEAlg

end # module
