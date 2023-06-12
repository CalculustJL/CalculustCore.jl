"""
$README
"""
module CalculustCore

using DocStringExtensions

include("Domains/Domains.jl")
include("Spaces/Spaces.jl")
include("BoundaryConditions/BoundaryConditions.jl")

using Adapt
using Functors
using Tricks
using SparseArrays

include("adapt.jl")
export cpu, gpu

using Reexport
@reexport using .Domains
@reexport using .Spaces
@reexport using .BoundaryConditions

#include("semidiscr.jl") # method of lines
#include("EigenValueProblem.jl")

const USE_CUDA = Ref{Union{Nothing, Bool}}(nothing)

@static if !isdefined(Base, :get_extension)
    import Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        Requires.@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
            include("../ext/CalculustCorePlotsExt.jl")
        end

        Requires.@require CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba" begin
            include("../ext/CalculustCoreCUDAExt.jl")
        end
    end
end

end # module
