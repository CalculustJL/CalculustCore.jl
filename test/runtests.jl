using AbstractPDEs
using Test, SafeTestsets

@testset "AbstractPDEs.jl" begin

    @time @safetestset "Field" begin include("Field.jl") end
    @time @safetestset "Operators" begin include("Operators.jl") end

#   @time @safetestset "Examples" begin include("examples.jl") end
end
