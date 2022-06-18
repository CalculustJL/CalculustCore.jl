#
using PDEInterfaces
using PDEInterfaces.SciMLOperators
using PDEInterfaces.LinearSolve

N = 32

@testset "1D Laplace" begin
    dom = reference_box(1)
    space = GaussLobattoLegendre1D(N; domain=dom)
    (x,) = pts = grid(space)

    @testset "Homogeneous Dirichlet-Dirichelt" begin
        op = laplaceOp(space)
        f  = @. 0*x + 1
        bcs = Dict(
                   :Lower1 => DirichletBC(),
                   :Upper1 => DirichletBC(),
                  )
        
        prob = BVPDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

        @time sol = solve(prob, alg; verbose=false)
        @test sol.resid < 1e-8
        plt = plot(sol)
        savefig(plt, "bvp_dd")
    end
    
    @testset "Homogeneous Dirichlet-Neumann" begin
        op = laplaceOp(space)
        f  = @. 0*x + 1
        bcs = Dict(
                   :Lower1 => DirichletBC(),
                   :Upper1 => NeumannBC(),
                  )
        
        prob = BVPDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

        @time sol = solve(prob, alg; verbose=false)
        @test sol.resid < 1e-8
        plt = plot(sol)
        savefig(plt, "bvp_dn")
    end
end

@testset "2D Laplace" begin
    dom = reference_box(2)
    space = GaussLobattoLegendre1D(N, N; domain=dom)
    (x, y,) = pts = grid(space)

    @testset "Homogeneous Dirichlet" begin
        op = laplaceOp(space)
        f  = @. sin(4π*x) * sin(3π*y)
#       f  = @. x * 0 + 1
        bcs = Dict(
                   :Lower1 => DirichletBC(),
                   :Upper1 => DirichletBC(),

                   :Lower2 => DirichletBC(),
                   :Upper2 => DirichletBC(),
                  )

        prob = BVPDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_GMRES())

        @time sol = solve(prob, alg; verbose=false)
        @test sol.resid < 1e-8
        plt = plot(sol)
        savefig(plt, "bvp2d_dd")
    end

    @testset "Dirichlet-Neumann" begin
        op = laplaceOp(space)
        f  = @. 0*x + 1
        bcs = Dict(
                   :Lower1 => NeumannBC(),
                   :Upper1 => DirichletBC(),

                   :Lower2 => DirichletBC(),
                   :Upper2 => NeumannBC(),
                  )

        prob = BVPDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

        @time sol = solve(prob, alg; verbose=false)
        @test sol.resid < 1e-8
        plt = plot(sol)
        savefig(plt, "bvp2d_dn")
    end
end
#
