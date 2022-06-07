#
using PDEInterfaces, LinearAlgebra
using PDEInterfaces.LinearSolve
using PDEInterfaces:solve
using PDEInterfaces.SciMLOperators

@testset "1D GLL" begin

    N = 128
    domain = reference_box(1)
    space = GaussLobattoLegendre1D(N; domain=domain)
    (x,) = grid = get_grid(space)

    D = gradOp(space) |> first
    M = massOp(space)
    A = laplaceOp(space)

    xp = range(-1,1;length=1000) |> Array
    J  = PDEInterfaces.lagrange_interp_mat(xp, x)

    u0 = @. 0*x + 1
    u1 = @. 1.0*x
    u2 = @. sin(10π * x)
    u4 = @. exp(x) / cos(x)
    u5 = @. 1 / (1+16x^2)

    # interpolation
    u4p = @. exp(xp) / cos(xp)
    u5p = @. 1 / (1+16xp^2)
    @test ≈(u4p, J * u4; atol=1e-8)
    @test ≈(u5p, J * u5; atol=1e-8)

    # Gradient
    @test ≈(D * u0, @. 0.0*x; atol=1e-8)
    @test ≈(D * u1, @. 0.0*x + 1.0; atol=1e-8)
    @test ≈(D * u2, @. (10π) * cos(10π * x); atol=1e-8)

    # Mass
    @test ≈(sum(M * u0), 2.0; atol=1e-8)
    @test ≈(sum(M * u1), 0.0; atol=1e-8)
    @test ≈(sum(M * u2), 0.0; atol=1e-8)

    # Laplacian
    Id = Diagonal([true for i=1:N])
    R  = Id[2:end-1, :] |> MatrixOperator
    AA = R * A * R'
    MM = R * M

    bb = MM * u2
#   sol = solve(LinearProblem(AA, bb), IterativeSolversJL_CG(); verbose=true)
    uu = AA \ bb
    u  = R' * uu
    @test ≈(-u, -(1 / (10π)^2) .* u2; atol=1e-8)
end
#
