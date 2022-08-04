
"""
1D Burgers + NN forcing

∂t(u) + u∂x(u) = νΔ(u) + NN(u)
"""
struct Model0 <: AbstractModel
    datafile
    data
    odeproblem
    space
    λ
end


function Model0(datafile::String)
end
