#
using PDEInterfaces, LinearSolve
using Plots

#=
@testset "1D Laplace" begin
    domain = reference_box(1)
    space = GaussLobattoLegendre1D(32; domain=domain)
    (x,) = grid = get_grid(space)

    @testset "Homogeneous Dirichlet-Dirichelt" begin
        op = laplaceOp(space)
        f  = @. 0*x + 1
        bcs = Dict(
                   :Lower1 => DirichletBC(),
                   :Upper1 => DirichletBC(),
                  )
        
        prob = BoundaryValuePDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_GMRES())
        u = solve(prob, alg)
        plt = plot(x, u)
        savefig(plt, "bvp_dd")
    end
    
    @testset "Homogeneous Dirichlet-Neumann" begin
        op = laplaceOp(space)
        f  = @. 0*x + 1
        bcs = Dict(
                   :Lower1 => DirichletBC(),
                   :Upper1 => NeumannBC(),
                  )
        
        prob = BoundaryValuePDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_GMRES())
    
        u = solve(prob, alg)
        plt = plot(x, u)
        savefig(plt, "bvp_dn")
    end
end
=#

@testset "2D Laplace" begin
    domain = reference_box(2)
    space = GaussLobattoLegendre1D(32, 32; domain=domain)
    (x, y,) = grid = get_grid(space)

    @testset "Homogeneous Dirichlet" begin
        op = laplaceOp(space)
        f  = @. 0*x + 1
        bcs = Dict(
                   :Lower1 => DirichletBC(),
                   :Upper1 => DirichletBC(),

                   :Lower2 => DirichletBC(),
                   :Upper2 => DirichletBC(),
                  )

        prob = BoundaryValuePDEProblem(op, f, bcs, space)
        alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_GMRES())
        u = solve(prob, alg)
        plt = plot(x, u)
        savefig(plt, "bvp2d_dd")
    end
end
#
