#
using Test, SafeTestsets

@testset "AbstractPDEInterfaces.jl" begin
    @safetestset "Domain" begin
        #@safetestset "name" begin include("file.jl") end
    end

    @safetestset "Spaces" begin
        #@safetestset "name" begin include("file.jl") end
    end

    @safetestset "Boundary Value Problem" begin
        #@safetestset "name" begin include("file.jl") end
    end

    #@time @safetestset "Examples" begin include("examples.jl") end
end
