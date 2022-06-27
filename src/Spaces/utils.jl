#
include("NDgrid.jl")
using SciMLOperators: ⊗, IdentityOperator, _reshape, _vec

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
