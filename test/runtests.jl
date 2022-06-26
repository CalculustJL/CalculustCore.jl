#
using Test, SafeTestsets

@testset "PDEInterfaces.jl" begin

@safetestset "Spaces" begin
#   @safetestset "Lagrange Polynomials" begin include("lagrangepoly.jl") end
    @safetestset "Fourier" begin include("fourier.jl") end
end

#@safetestset "Boundary Value Problem" begin
#    @safetestset "BVP" begin include("bvp.jl") end
#end

@time @safetestset "Examples" begin include("examples.jl") end
end
