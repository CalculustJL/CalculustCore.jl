#
using PDEInterfaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(PDEInterfaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using JLD2, Plots

function process(name)
    data = jldopen(name)

    x = data["x"]
    t = data["t"]
    u = data["pred"]

    nx, ns, nt = size(u)

    for i=1:20:ns
        @views ut = u[:,i,:]
        anim8_traj(ut, t, x; name="traj"*string(i))
    end
end

function plot_traj(u, t, x)
    plt = plot()
    for i=1:length(t)
        plot!(plt, x, u[:,i], legend=false)
    end
    plt
end

function anim8_traj(u, t, x; name=nothing)
    ylims = begin
        u0 = @views u[:,1]
        mi = minimum(u0)
        ma = maximum(u0)
        buf = (ma-mi)/5
        (mi-buf, ma+buf)
    end

    anim = @animate for i=1:length(t)
        plt = plot(x, u[:,i], legend=false, ylims=ylims)
    end

    filename = joinpath(dirname(@__FILE__), name * ".gif")
    gif(anim, filename, fps=10)
end

#================#
name = "burgers_nu1em5_n8192"
filename = joinpath(@__DIR__, name * ".jld2") 

#process(filename)
nothing
#
