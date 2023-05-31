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

@testset "Unit Domains" begin
    for dom in (
        UnitIntervalDomain(),
        UnitSquareDomain(),
        UnitCubeDomain(),
        UnitBoxDomain(4),
    )

        D = dims(dom)
        @test isperiodic(dom) === Tuple(false for _ in 1:D)
        @test [expanse(dom)...] ≈ ones(D)

        bdr = boundaries(dom)
        @test length(bdr) == 2D

        for i in eachindex(bdr)
            _bdr = bdr[i]
            @test dims(_bdr) == D - 1
        end

        # TODO to test boundary locations, overload Base.in
        # TODO: have a notion of what space each domain is embedded in
        # eg {3} × (0, 1) ∈ R^3

        @test domain_tag(dom) === :Interior

        if D == 1
            @test boundary_tags(dom) === (:D1_inf, :D1_sup,)
        elseif D == 2
            @test boundary_tags(dom) === (:D1_inf, :D1_sup,
                                          :D2_inf, :D2_sup,)
        elseif D == 3
            @test boundary_tags(dom) === (:D1_inf, :D1_sup,
                                          :D2_inf, :D2_sup,
                                          :D3_inf, :D3_sup,)
        elseif D == 4
            @test boundary_tags(dom) === (:D1_inf, :D1_sup,
                                          :D2_inf, :D2_sup,
                                          :D3_inf, :D3_sup,
                                          :D4_inf, :D4_sup)
        end
    end
end

@testset "Fourier Domain" begin
    for D in 1:4

        dom = FourierDomain(D)

        @test isperiodic(dom) === Tuple(true for _ in 1:D)
        @test [expanse(dom)...] ≈ 2pi * ones(D)

        bdr = boundaries(dom)
        @test length(bdr) == 2D

        @test domain_tag(dom) === :Interior

        if D == 1
            @test boundary_tags(dom) === (:D1_periodic, :D1_periodic,)
        elseif D == 2
            @test boundary_tags(dom) === (:D1_periodic, :D1_periodic,
                                          :D2_periodic, :D2_periodic,)
        elseif D == 3
            @test boundary_tags(dom) === (:D1_periodic, :D1_periodic,
                                          :D2_periodic, :D2_periodic,
                                          :D3_periodic, :D3_periodic,)
        elseif D == 4
            @test boundary_tags(dom) === (:D1_periodic, :D1_periodic,
                                          :D2_periodic, :D2_periodic,
                                          :D3_periodic, :D3_periodic,
                                          :D4_periodic, :D4_periodic)
        end
    end
end

@testset "BoxDomain" begin
    x0, x1 = 2.0, 3.0
    y0, y1 = 5.0, 7.0
    z0, z1 = 4.0, 8.0

    bdr_tags = (:x0, :x1, :y0, :y1, :z0, :z1)

    dom = BoxDomain(x0, x1, y0, y1, z0, z1;
                    periodic_dims = (2, 3,),
                    tag = :xyz,
                    boundary_tags = bdr_tags)

    D = 3
    @test dims(dom) == D
    @test isperiodic(dom) === (false, true, true,)
    @test [expanse(dom)...] ≈ [x1-x0, y1-y0, z1-z0]

    bdr = boundaries(dom)
    @test length(bdr) == 2D

    for i in eachindex(bdr)
        _bdr = bdr[i]
        @test dims(_bdr) == D - 1
    end

    @test domain_tag(dom) === :xyz
    @test boundary_tags(dom) === bdr_tags
end
#
