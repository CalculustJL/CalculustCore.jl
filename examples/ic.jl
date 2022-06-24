
using Random

Random.seed!(0)

function uIC(x, ftr, k)
    u0 = @. sin(2x) + sin(3x) + sin(5x)
#   u0 = @. sin(x - π)

    u0 = begin
        u  = rand(size(x)...)
        uh = ftr * u

#       uh = rand(ComplexF64, size(k)...)
#       uh[1] = rand()
        uh[20:end] .= 0
        ftr \ uh
    end

    u0
end


plot(x, uIC(x, ftr, k))
#
