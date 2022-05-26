#
using PDEInterfaces
using Test, SafeTestsets

@testset "PDEInterfaces.jl" begin

    @time @safetestset "Field" begin include("Field.jl") end
    @time @safetestset "Operators" begin include("Operators.jl") end
    @time @safetestset "Boundary Value Problem" begin include("bvp.jl") end

#   @time @safetestset "Examples" begin include("examples.jl") end
end
