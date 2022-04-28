#
using AbstractPDEs, LinearAlgebra

nr = 8
ns = 12

arr = rand(nr,ns)
u   = Field(arr)

@test eachindex(vec(arr)) == eachindex(u)

@inferred similar(u)
@inferred +u
@inferred -u
@inferred 2u
@inferred 2 .+ u
@inferred 2 .- u
@inferred u .+ 2
@inferred u .- 2

