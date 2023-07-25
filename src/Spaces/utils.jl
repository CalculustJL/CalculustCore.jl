#
function _pair_update_funcs(vecs, funcs)
    VV = AbstractSciMLOperator[]

    for i in 1:length(vecs)
        func = if funcs isa Nothing
            DEFAULT_UPDATE_FUNC
        else
            funcs[i]
        end

        V = DiagonalOperator(vecs[i]; update_func = func)
        push!(VV, V)
    end

    VV
end

"""
(f1 ∘ f2)(v, u, p, t) with caching
"""
struct ComposedUpdateFunction{F1, F2, C <: AbstractArray}
    f1::F1
    f2::F2
    cache::C

    function ComposedUpdateFunction(f1, f2, cache)
        f1 = isnothing(f1) ? DEFAULT_UPDATE_FUNC : f1
        f2 = isnothing(f2) ? DEFAULT_UPDATE_FUNC : f2

        new{typeof(f1), typeof(f2), typeof(cache)}(f1, f2, cache)
    end
end

function (A::ComposedUpdateFunction)(u, p, t)

    v = A.f2(u, p, t)
    A.f1(v, p, t)
end

function (A::ComposedUpdateFunction)(v, u, p, t)

    A.f2(A.cache, u, p, t)
    A.f1(v, A.cache, p, t)
end
#
