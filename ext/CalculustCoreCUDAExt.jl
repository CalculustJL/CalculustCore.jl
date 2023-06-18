module CalculustCoreCUDAExt

using CUDA
using Adapt
using CalculustCore

function CalculustCore.check_use_cuda()
    CalculustCore.USE_CUDA[] = CUDA.functional()
end

## CUDA Adaptor
function Adapt.adapt_storage(::CalculustCore.CalculustCUDAAdaptor, x)
    cu(x)
end

## CPU Adaptor
function Adapt.adapt_storage(::CalculustCore.CalculustCPUAdaptor, x::CUDA.CUSPARSE.AbstractCuSparseMatrix,)
    adapt(Array, x)
end

end
