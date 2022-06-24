#
dir = "../examples"
files = (
         "lagrangepoly/1D_Poisson.jl",

         "fourier/advect.jl",
         "fourier/heat.jl",
         "fourier/heat_forcing.jl",
        )
for file in files
    @time @safetestset "$file" begin include("$dir/$file") end
end
#
