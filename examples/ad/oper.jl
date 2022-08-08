#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using LinearAlgebra, SciMLSensitivity, Zygote

N = 16

u0 = rand(N)
ps = rand(N)

space = FourierSpace(N)
space = make_transform(space, u0)

F  = transformOp(space)
Dx = gradientOp(space)[1]

S  = rand() * MatrixOperator(rand(N,N))
Di = DiagonalOperator(rand(N))
M  = MatrixOperator(rand(N,N))
Af = AffineOperator(rand(N,N), rand(N,N), rand(N))
Ad = MatrixOperator(rand(N,N)) + MatrixOperator(rand(N,N))
Co = MatrixOperator(rand(N,N)) * MatrixOperator(rand(N,N))
T  = TensorProductOperator(rand(n,n), rand(n,n))
Id = IdentityOperator{N}()
Z  = NullOperator{N}()
Ao = SciMLOperators.AdjointOperator(rand(N,N) |> MatrixOperator)
Tr = SciMLOperators.TransposedOperator(rand(N,N) |> MatrixOperator)

α  = ScalarOperator(rand())
β  = ScalarOperator(rand()) + ScalarOperator(rand())
γ  = ScalarOperator(rand()) * ScalarOperator(rand())

# multidim

loss = function(p)

    v = Diagonal(p) * u0
    #v = Zygote.hook(Δ -> (println("Δv: ", Δ); Δ), v)
    v = Zygote.hook(Δ -> (println("Δv: ", typeof(Δ)); Δ), v)

    #w = Dx * v ## Δ vanishes - ComposedOperator # INCORRECT
    w = F \ F * v ## Δ vanishes - ComposedOperator # INCORRECT
    #w = Co * v ## Δ vanishes - ComposedOperator # INCORRECT
    #w = β  * v ## Δ - AddedScalarOperator # ERROR

    #w = S * v   ## Δ ok - ScaledOperator
    #w = T * v   ## Δ ok - TensorProductOperator
    #w = Ad * v  ## Δ ok - AddedOperator
    #w = Af * v  ## Δ ok - AffineOperator
    #w = Di * v  ## Δ ok - DiagonalOperator
    #w = M  * v  ## Δ ok - MatrixOperator
    #w = Id * v  ## Δ ok - IdentityOperator
    #w = Z  * v  ## Δ ok - NullOperator
    #w = Ao' * v  ## Δ ok - AdjointOperator
    #w = transpose(Tr) * v  ## Δ ok - TransposedOperator
    #w = F\(F*v) ## Δ ok - FunctionOperator

    #w = α  * v  ## Δ ok - ScalarOperator
    #w = γ  * v  ## Δ ok - ComposedScalarOperator
    w = Zygote.hook(Δ -> (println("Δw: ", Δ); Δ), w)

    l = sum(w)
    l = Zygote.hook(Δ -> (println("Δl: ", Δ); Δ), l)
end

# dummy calls
println("fwd"); @time loss(ps) |> display                  # works
println("bwd"); @time Zygote.gradient(loss, ps) |> display # works
#
