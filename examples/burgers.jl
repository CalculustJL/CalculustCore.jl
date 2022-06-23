#
# add dependencies to env stack
pkgpath = dirname(dirname(@__FILE__))
tstpath = joinpath(pkgpath, "test")
!(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)

using PDEInterfaces
using OrdinaryDiffEq, LinearSolve, Sundials
using Plots

N = 4096
ν = 1e-3
p = ()

function uIC(x, ftr)
    u0 = @. sin(2x) + sin(3x) + sin(5x)
#   u0 = @. sin(x - π)
    
#   u0 = begin
#       uh = rand(ComplexF64, size(k))
#       uh[20:end] .= 0
#       ftr \ uh
#   end

    u0
end
    
function solve_burgers(N, ν, p; uIC=uIC)

    """ space discr """
    space = FourierSpace(N)
    discr = Collocation()
    
    (x,) = points(space)
    ftr  = transforms(space)
    k = modes(space)

    """ IC """
    u0 = uIC(x, ftr)
    
    """ operators """
    A = diffusionOp(ν, space, discr)
    
    function burgers!(v, u, p, t)
        copy!(v, u)
        v
    end
    
    v = @. x*0 + 1
    f = @. x*0
    C = advectionOp((v,), space, discr; vel_update_funcs=(burgers!,))
    
    F = AffineOperator(-C, f)

    A = cache_operator(A, x)
    F = cache_operator(F, x)
    
    """ time discr """
    tspan = (0.0, π)
    tsave = range(tspan...; length=10)
    
    #odealg = Rodas5(autodiff=false)
    #odealg = Tsit5()
    odealg = CVODE_BDF(method=:Functional)
    
    prob = SplitODEProblem(A, F, u0, tspan, p)
    @time sol = solve(prob, odealg, saveat=tsave)
    
    sol
end

sol = solve_burgers(N, ν, p)

""" analysis """
plt = plot()
for i=1:length(sol.u)
    plot!(plt, x, sol.u[i], legend=false)
end
display(plt)

nothing
