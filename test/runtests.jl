#
using Test, SafeTestsets

@testset "CalculustCore.jl" begin
    @safetestset "Domain" begin
        @safetestset "name" begin include("domains.jl") end
    end

    @safetestset "Spaces" begin
        #@safetestset "name" begin include("file.jl") end
    end

    @safetestset "Boundary Value Problem" begin
        #@safetestset "name" begin include("file.jl") end
    end
end
