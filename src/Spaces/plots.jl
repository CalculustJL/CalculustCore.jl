#
"""
plot of function over space
args:
    - u scalar field
    - space AbstractSpace
"""
function Plots.plot(u::AbstractVector, space::AbstractSpace{<:Any, 1}; kwargs...)
    (x,) = points(space)
    plt = plot(x, u; kwargs...)
end

function Plots.plot(u::AbstractVector, space::AbstractSpace{<:Any, 2}; a = 30, b = 30,
                    kwargs...)
    npts = size(space)
    (x, y) = points(space)

    u = reshape(u, npts)
    x = reshape(x, npts)
    y = reshape(y, npts)

    #   plt = plot(x, y, u, legend=false, c=:grays, camera=(a,b))
    #   plt = plot!(x', y', u', legend=false, c=:grays, camera=(a,b))

    plt = Plots.heatmap(u; kwargs...)

    plt
end

function Plots.animate(u::AbstractMatrix,
                       space::AbstractSpace{<:Any, 1},
                       t::AbstractVector = collect(1:size(u, 2));
                       kwargs...)
    ylims = begin
        mi = minimum(u)
        ma = maximum(u)
        buf = (ma - mi) / 5
        (mi - buf, ma + buf)
    end
    anim = @animate for i in 1:size(u, 2)
        titlestr = "time = $(round(t[i], digits=8))"
        plt = plot(u[:, i], space; ylims = ylims, title = titlestr, kwargs...)
    end
end

function Plots.animate(u::AbstractMatrix,
                       space::AbstractSpace{<:Any, 2},
                       t::AbstractVector = collect(1:size(u, 2));
                       kwargs...)
    clim = begin
        mi = minimum(u)
        ma = maximum(u)
        (mi, ma)
    end
    anim = @animate for i in 1:size(u, 2)
        titlestr = "time = $(round(t[i], digits=8))"
        plt = plot(u[:, i], space; clim = clim, title = titlestr, kwargs...)
    end
end
