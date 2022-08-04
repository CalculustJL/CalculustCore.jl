
module Burgers1D

using PDEInterfaces

using OrdinaryDiffEq, CUDA, LinearAlgebra, ComponentArrays
using Lux, Random, JLD2, SciMLSensitivity, Zygote
using Optimization, OptimizationOptimJL, OptimizationOptimisers, Optimisers
using Plots

abstract type AbstractModel end

function predict(::AbstractModel) end
function loss(::AbstractModel) end
function save(::AbstractModel) end

include("Model0.jl")
include("Model1.jl")

###
# utils
###

# get stuff from data
function ut_from_data(datafile)
    data = jldopen(datafile)
    
    t = data["t"]
    u = data["u_coarse"]

    u, t
end

# callbacks
function optcb(p, l, pred;
                doplot=false,
                space=space,
                steptime=nothing,
                iter=nothing,
                niter=nothing,
               )

    steptime = steptime isa Nothing ? 0.0 : steptime
    iter = iter isa Nothing ? 0 : iter
    niter = niter isa Nothing ? 0 : niter

    println(
            "[$iter/$niter] \t Time $(round(steptime; digits=2))s \t Loss: " *
            "$(round(l; digits=8)) \t "
           )

    return false
end

# training
function train(loss, p;
               alg=Optimisers.Adam(1f-2),
               maxiters=1000,
               callback=optcb,
              )

    adtype = Optimization.AutoZygote()
    # x=object to optimize
    # p=parameters for optimization loop
    optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, p)

    optres = Optimization.solve(optprob, alg, callback=callback, maxiters=maxiters)

    optres.u
end

end #module
