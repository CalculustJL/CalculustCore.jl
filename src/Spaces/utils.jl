#
include("NDgrid.jl")
import SciMLOperators: _reshape, _vec, copy!

_reshape(a::CUDA.CuArray, dims::NTuple{D,Int}) where{D} = reshape(a,dims)

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

###
# GPU
###

import Lux: cpu, gpu, LuxCPUAdaptor, LuxCUDAAdaptor

_fft_lib(u::AbstractArray) = FFTW
_fft_lib(u::CUDA.CuArray) = CUDA.CUFFT


##
## http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga9f64da9e0125c712672fd89d166d3b9c.html
#import LinearAlgebra: mul!
#function LinearAlgebra.mul!(v::CuVector, D::Diagonal{<:Any,CuVector}, u::CuVector)
#    _copy!(v, u)
#    n = length(D)
#    CUDA.CUBLAS.tbmv!(
#                      'u', # upper triangular
#                      'n', # transpose false
#                      'n', # not unit triangular
#                      0,   # k - num super diagonals
#                      D.diag,
#                      v
#                     )
#end
#
