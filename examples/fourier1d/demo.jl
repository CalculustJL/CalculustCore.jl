using PDEInterfaces, OrdinaryDiffEq, Plots
N = 1024; ν = 1e-3; tspan=(0.0,10.0)

""" space discr """
space = FourierSpace(N)
(x,) = points(space)

burgers!(v, u, p, t) = copy!(v, u)
forcing!(f, u, p, t) = lmul!(false, f)

A = diffusionOp(ν, space, Collocation())
C = advectionOp((zero(x),), space, Collocation(); vel_update_funcs=(burgers!,))
F = -C + forcingOp(zero(x), space, Collocation(); f_update_func=forcing!)

A = cache_operator(A, x)
F = cache_operator(F, x)

u0 = truncationOp(space; truncation_frac=1//4) * rand(N)
prob = SplitODEProblem(A, F, u0, tspan, p; reltol=1e-8, abstol=1e-8)
@time sol = solve(prob, odealg, saveat=range(tspan...;length=100))

anim = animate(Array(sol), space)
gif(anim, "burgers.gif", fps= 20)






