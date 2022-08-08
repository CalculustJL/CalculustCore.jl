#
include("NDgrid.jl")

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

## update function composition
"""
f1 ∘ f2
"""
struct ComposedUpdateFunction{F1,F2}
    f1::F1
    f2::F2

    function ComposedUpdateFunction(f1 = nothing, f2 = nothing)
        id = function(v, u, p, t)
            v
        end

        f1 = f1 isa Nothing ? id : f1
        f2 = f2 isa Nothing ? id : f2

        new{typeof(f1),typeof(f2)}(f1, f2)
    end
end

function (C::ComposedUpdateFunction)(v, u, p, t)
    C.f2(v, u, p, t)
    C.f1(v, u, p, t)
end

###
# GPU
###

import Lux: cpu, gpu, LuxCPUAdaptor, LuxCUDAAdaptor

_fft_lib(u::AbstractArray) = FFTW
_fft_lib(u::CUDA.CuArray) = CUDA.CUFFT
#
