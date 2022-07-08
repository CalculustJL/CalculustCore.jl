#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using OrdinaryDiffEq, LinearSolve, LinearAlgebra, Random
using Plots

nx = 128
ny = 128
ν = 1e-2
p = nothing

Random.seed!(0)
function uIC(space)
    x = points(space)[1]
    X = truncationOp(space,1//20)

    u0 = X * rand(size(x)...)

    u0
end

function solve_burgers(nx, ny, ν, p;
                       uIC=uIC,
                       tspan=(0.0, 10.0),
                       nsave=100,
                       odealg=SSPRK43(),
                      )

    """ space discr """
    space = FourierSpace(nx, ny)
    discr = Collocation()

    x, y = points(space)

    """ IC """
    u0 = uIC(space)

    """ operators """
    A = diffusionOp(ν, space, discr)

    function burgers!(v, u, p, t)
        copy!(v, u)
    end

    function forcing!(f, u, p, t)
        lmul!(false, f)
    end

    C = advectionOp((zero(x),), space, discr; vel_update_funcs=(burgers!,))
    F = -C + forcingOp(zero(x), space, discr; f_update_func=forcing!)

    A = cache_operator(A, x)
    F = cache_operator(F, x)

    """ time discr """
    function Ajac(Jv, v, u, p, t;A=A)
        SciMLOperators.update_coefficients!(A, u, p, t)
        mul!(Jv, A, v)
    end
    odefunc = SplitFunction(A, F; jvp=Ajac)

    tsave = range(tspan...; length=nsave)
    prob = ODEProblem(odefunc, u0, tspan, p; reltol=1e-8, abstol=1e-8)
    @time sol = solve(prob, odealg, saveat=tsave)

    sol, space
end

sol, space = solve_burgers(nx, ny, ν, p)
pred = Array(sol)
anim = animate(pred, space)
gif(anim, "examples/fourier/a.gif", fps= 20)
#
