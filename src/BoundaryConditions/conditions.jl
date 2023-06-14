#
###
# Boundary Condition Types
###

"""
$TYPEDEF

u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct DirichletBC{F}
    f::F = (points...) -> zero(first(points))
end

"""
$SIGNATURES

(n⋅∇)u(x) = f(x), x ∈ ∂Ω

Defaults to homogeneous.
"""
Base.@kwdef struct NeumannBC{F}
    f::F = (points...) -> zero(first(points))
end

"""
$SIGNATURES

f1(x)u(x) + f2(x)(n⋅∇)u(x) = f3(x), x ∈ ∂Ω
"""
struct RobinBC{F1, F2, F3}
    f1::F1
    f2::F2
    f3::F3
end

"""
$SIGNATURES

Periodic Boundary Condition
"""
Base.@kwdef struct PeriodicBC{Ttag}
    tag::Ttag
end
#
