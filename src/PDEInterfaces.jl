module PDEInterfaces

using Reexport

#@reexport using SciMLBase
import SciMLBase; import SciMLBase:solve, init; export solve, init

using LinearSolve

import Plots: plot, plot!
import UnPack: @unpack

include("Domains/Domains.jl")
@reexport using Domains
using Domains: AbstractDomain

include("Spaces/Spaces.jl")
@reexport using Spaces
using Space: AbstractSpace

include("utils.jl")
include("BoundaryConditions.jl")
#include("Discretizations.jl")

# problems
include("BoundaryValueProblem.jl")
#include("ODEProblem.jl")
#include("EigenValueProblem.jl")

export 
       ## Fields
#      Field,

       ## Boundary conditions
       DirichletBC, NeumannBC, RobinBC, PeriodicBC,

       ## Boundary vale problem
       BVPDEProblem, LinearBVPDEAlg, NonlinearBVPDEAlg

end # module
