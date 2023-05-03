module AbstractPDEInterfacesCUDAExt

if isdefined(Base, :get_extension)
    using CUDA
    using AbstractPDEInterfaces
else
    using ..CUDA
    using ..AbstractPDEInterfaces
end

AbstractPDEInterfaces.USE_CUDA[] = CUDA.functional()

## CUDA Adaptor
function Adapt.adapt_storage(::AbstractPDEInterfaces.CUDAAdaptor, x)
    Adapt.adapt_storage(CUDA.Adaptor, x) # cu(x)
end

## CPU Adaptor
function Adapt.adapt_storage(::CPUAdaptor, x::CUDA.CUSPARSE.AbstractCuSparseMatrix,)
    adapt(Array, x)
end

end
