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

    pred  = data["pred"]
    time  = data["time"]
    space = data["space"]

    npts, nsims, ntime = size(pred)
    x = points(space) |> first

    for i=1:5:nsims
        @views upred = pred[:,1,:]
    end

end

function plot_traj(u, t, x)
    plt = plot()
    for i=1:length(time)
        plot!(plt, x, pred[:,i], legend=false)
    end
    plt
end

function anim8_traj(pred, time, x)
    ylims = begin
        u = @views pred[:,1]
        mi = minimum(u)
        ma = maximum(u)
        buf = (ma-mi)/5
        (mi-buf, ma+buf)
    end
    anim = @animate for i=1:length(time)
        plt = plot(x, pred[:,i], legend=false, ylims=ylims)
    end
end

#
