#
using Zygote, LinearAlgebra

N = 4
u0 = rand(N)
ps = rand(N)

mats = (rand(N,N), rand(N,N),) # (A, B,)
nums = (rand(), rand(),)       # (α, β,)

loss_m = function(p)
    v = Diagonal(p) * u0
    v = Zygote.hook(Δ -> (println("Δv: ", typeof(Δ)); Δ), v)

    w = foldl((acc, op) -> op * acc, mats; init=v) # w = B * A * v
    w = Zygote.hook(Δ -> (println("Δw: ", Δ); Δ), w)

    l = sum(w)
    l = Zygote.hook(Δ -> (println("Δl: ", Δ); Δ), l)
end

println("fwd"); @time loss_m(ps) |> display
println("bwd"); @time Zygote.gradient(loss_m, ps) |> display # INCORRECT - should not vanish

loss_n = function(p)
    v = Diagonal(p) * u0
    v = Zygote.hook(Δ -> (println("Δv: ", typeof(Δ)); Δ), v)

    w = sum(a -> convert(Number, a), nums; init=zero(eltype(nums))) * v # w = αβ * v
    #w = sum(a -> convert(Number, a), nums) * v # w = αβ * v
    w = Zygote.hook(Δ -> (println("Δw: ", Δ); Δ), w)

    l = sum(w)
    l = Zygote.hook(Δ -> (println("Δl: ", Δ); Δ), l)
end

println("fwd"); @time loss_n(ps) |> display
println("bwd"); @time Zygote.gradient(loss_n, ps) |> display # ERRORS

#
