#
using Test
using CalculustCore
using CalculustCore.Domains: NullDomain, ∅, PointDomain, ProductDomain

# NullDomain
@test_throws ArgumentError bounding_box(∅)
@test_throws ArgumentError expanse(∅)
@test_throws ArgumentError isperiodic(∅)
@test boundaries(∅) === ()
@test_throws ArgumentError domain_tag(∅)

# PointDomain

for dom in (
            UnitIntervalDomain(),
            UnitSquareDomain(),
            UnitCubeDomain(),
           )

    D = dims(dom)
    @test isperiodic(dom) === Tuple(false for _ in 1:D)
    @test [expanse(dom)...] ≈ ones(D)

    bdr = boundaries(dom)
    @test length(bdr) == 2D

    # test boundary locations

    @test domain_tag(dom) === :Interior

    if D == 1
        @test boundary_tags(dom) === (:D1_inf, :D1_sup,)
    elseif D == 2
        @test boundary_tags(dom) === (:D1_inf, :D1_sup, :D2_inf, :D2_sup,)
    elseif D == 3
        @test boundary_tags(dom) === (:D1_inf, :D1_sup, :D2_inf, :D2_sup, :D3_inf, :D3_sup,)
    end
end

