#
abstract type DeviceAdaptor end
struct CPUAdaptor <: DeviceAdaptor end
struct CUDAAdaptor <: DeviceAdaptor end

Adapt.adapt_storage(::CPUAdaptor, x::AbstractArray) = adapt(Array, x)
function Adapt.adapt_storage(::CPUAdaptor, x::Union{AbstractRange, SparseArrays.AbstractSparseArray},)
    x
end

_isbitsarray(::AbstractArray{<:Number}) = true
_isbitsarray(::AbstractArray{T}) where{T} = isbitstype(T)
_isbitsarray(x) = false

_isleaf(x) = _isbitsarray(x) || Functors.isleaf(x)

"""
    cpu(x)

Transfer `x` to CPU
"""
function cpu end

cpu(x) = fmap(x -> adapt(CPUAdaptor(), x), x)

"""
    gpu(x)

Transfer `x` to GPU
"""
function gpu(x)
    has_cuda = static_hasmethod(check_use_cuda, typeof(()))

    if !has_cuda
        @warn "CUDA is not loaded. Defaulting to CPU."
        return x
    end

    check_use_cuda()
    if USE_CUDA[]
        # return adapt(CUDAAdaptor(), x)
        return fmap(x -> adapt(CUDAAdaptor(), x), x; exclude=_isleaf)
    else
        @warn "CUDA loaded but not CUDA.functional() = false. Defaulting to CPU."
        return x
    end
end

function check_use_cuda end
