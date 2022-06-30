#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using CUDA, LinearAlgebra, BenchmarkTools
using OrdinaryDiffEq

p = nothing
t = 0f0

K = 10
for N in 2 .^(14:14)
    println("problem size, N: ", N)
    println("num cases, K: ", K)

    space = FourierSpace(N) |> gpu
    (x,)  = points(space)
    discr = Collocation()

    function burgers!(v, u, p, t)
        copyto!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
    end

    A = diffusionOp(1f0, space, discr)

    C = advectionOp((zero(x),), space, discr; vel_update_funcs=(burgers!,))
    F = -C + forcingOp(zero(x), space, discr; f_update_func=forcing!)

#   odefunc = SplitFunction(A, F)

    println("########################")
    println("Baseline, [u]")
    println("########################")

    let
        u  = copy(x)
        du = copy(x)

        A = cache_operator(A, u)
        F = cache_operator(F, u)

        println("imp iip")
        @btime $A($du, $u, $p, $t)

#       println("imp oop")
#       @btime $A($u, $p, $t)

        println("exp iip")
        @btime $F($du, $u, $p, $t)

#       println("exp oop")
#       @btime $F($u, $p, $t)
    end

    println("########################")
    println("Batching via Broadcast [u... K times])") # linear scaling
    println("########################")

    let
        u  = [copy(x) for i=1:K]
        du = [copy(x) for i=1:K]

        A = cache_operator(A, u[1])
        F = cache_operator(F, u[1])

        println("imp iip")
        @btime $A.($du, $u, $p, $t)

#       println("imp oop")
#       @btime $A.($u, $p, $t)

        println("exp iip")
        @btime $F.($du, $u, $p, $t)

#       println("exp oop")
#       @btime $F.($u, $p, $t)
    end

    println("########################")
    println("Batching via mul (N,K) array") # insane!!
    println("########################")

    let
        o  = ones(1,K) |> gpu
        u  = copy(x) * o
        du = copy(x) * o

        space = make_transform(space, u)

        A = diffusionOp(1f0, space, discr)

        C = advectionOp((zero(u),), space, discr; vel_update_funcs=(burgers!,))
        F = -C + forcingOp(zero(u), space, discr; f_update_func=forcing!)

        A = cache_operator(A, u)
        F = cache_operator(F, u)

        println("imp iip")
        @btime $A($du, $u, $p, $t)

#       println("imp oop")
#       @btime $A($u, $p, $t)

        println("exp iip")
        @btime $F($du, $u, $p, $t)

#       println("exp oop")
#       @btime $F($u, $p, $t)
    end

end
#
