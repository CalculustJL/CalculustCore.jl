module AbstractPDEInterfaces

using DocStringExtensions

include("Domains/Domains.jl")
include("Spaces/Spaces.jl")
include("BoundaryConditions/BoundaryConditions.jl")

using Reexport
@reexport using .Domains
@reexport using .Spaces
@reexport using .BoundaryConditions

#include("semidiscr.jl") # method of lines
#include("EigenValueProblem.jl")

end # module
