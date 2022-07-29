#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearAlgebra, Random
using Plots

N = 1024
ν = 1e-3
p = nothing

Random.seed!(0)
function uIC(space)
    x = points(space)[1]
    X = truncationOp(space, (1/8,))

    u0 = X * rand(size(x)...)
end

odecb = begin
    function affect!(int)
        println(
                "[$(int.iter)] \t Time $(round(int.t; digits=8))s"
               )
    end

    DiscreteCallback((u,t,int) -> true, affect!, save_positions=(false,false))
end

function solve_burgers1D(N, ν, p;
                         uIC=uIC,
                         tspan=(0.0, 10.0),
                         nsims=1,
                         nsave=100,
                         odealg=SSPRK43(),
                        )

    """ space discr """
    space = FourierSpace(N)
    discr = Collocation()

    (x,) = points(space)

    """ IC """
    u0 = [uIC(space) for i=1:nsims]
    u0 = hcat(u0...)
    space = make_transform(space, u0; p=p)

    """ operators """
    function burgers!(v, u, p, t)
        copyto!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
    end

    A = diffusionOp(ν, space, discr)
    C = advectionOp((zero(u0),), space, discr; vel_update_funcs=(burgers!,))
    F = forcingOp(zero(u0), space, discr; f_update_func=forcing!)
    odefunc = cache_operator(A-C+F, u0) |> ODEFunction

    """ time discr """
    tsave = range(tspan...; length=nsave)
    prob = ODEProblem(odefunc, u0, tspan, p; reltol=1e-8, abstol=1e-8)
    @time sol = solve(prob, odealg, saveat=tsave, callback=odecb)

    sol, space
end

sol, space = solve_burgers1D(N, ν, p)
space = cpu(space)
pred = Array(sol)
nothing
#
