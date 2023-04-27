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

import Requires

@static if !isdefined(Base, :get_extension)
    function __init__()
        Requires.@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
            include("../ext/AbstractPDEInterfacesPlotsExt.jl")
        end
    end
end

end # module
