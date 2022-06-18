#
dir = "../examples"
files = (
         "p2d.jl",
        )
for file in files
    @time @testset "$file" begin include("$dir/$file") end # @safetestset ?
end
#
