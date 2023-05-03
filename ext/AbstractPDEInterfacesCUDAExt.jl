module AbstractPDEInterfacesCUDAExt

if isdefined(Base, :get_extension)
    using CUDA
    using Adapt
    using AbstractPDEInterfaces
else
    using ..CUDA
    using ..Adapt
    using ..AbstractPDEInterfaces
end

function AbstractPDEInterfaces.check_use_cuda()
    AbstractPDEInterfaces.USE_CUDA[] = CUDA.functional()
end

## CUDA Adaptor
function Adapt.adapt_storage(::AbstractPDEInterfaces.CUDAAdaptor, x)
    cu(x)
end

## CPU Adaptor
function Adapt.adapt_storage(::AbstractPDEInterfaces.CPUAdaptor, x::CUDA.CUSPARSE.AbstractCuSparseMatrix,)
    adapt(Array, x)
end

end
