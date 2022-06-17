#
using Test, SafeTestsets

@testset "PDEInterfaces.jl" begin

@time @safetestset "Vector Calculus Operators" begin include("vec_calc.jl") end
#@time @safetestset "Boundary Value Problem" begin include("bvp.jl") end
#@time @safetestset "Examples" begin include("examples.jl") end
end
