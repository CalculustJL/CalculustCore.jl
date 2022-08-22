#
_transp(a, ::AbstractDiscretization) = transpose(a)

function _pair_update_funcs(vecs, funcs)

    VV = AbstractSciMLOperator[]

    for i=1:length(vecs)
        func = if funcs isa Nothing
            DEFAULT_UPDATE_FUNC
        else
            funcs[i]
        end

        V = DiagonalOperator(vecs[i]; update_func=func)
        push!(VV, V)
    end

    VV
end

"""
(f1 ∘ f2)(v, u, p, t) with caching
"""
struct ComposedUpdateFunction{F1,F2,C<:AbstractArray}
    f1::F1
    f2::F2
    cache::C

    function ComposedUpdateFunction(f1 = nothing, f2 = nothing, cache=nothing)
        f1 = f1 isa Nothing ? DEFAULT_UPDATE_FUNC : f1
        f2 = f2 isa Nothing ? DEFAULT_UPDATE_FUNC : f2

        new{typeof(f1),typeof(f2), typeof(cache)}(f1, f2, cache)
    end
end

function (A::ComposedUpdateFunction)(v, u, p, t)
    @unpack f1, f2, cache = A

    f2(cache, u, p, t)
    f1(v, cache, p, t)
end

#
