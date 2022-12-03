ESI = ExactSparseInterpolations

maxdeg(f) = maximum(map(t -> total_degree(t), terms(f)), init=0) 

@testset "van-der-Hoeven-Lecerf, Q" begin
    R, (x1,x2) = PolynomialRing(QQ, ["x1","x2"])
    cases = [
        (x1)//(2x2 + 1), R(2)//R(3), R(0)//R(1),
        (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3,
        (x1 + x2)^5//(x1 - x2)^6, (x1*x2 - x1 - x2)^3 // R(8)
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        vdhl = ESI.vanDerHoevenLecerf(R, n, d)
        P, Q = ESI.interpolate!(vdhl, f)
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

    f = ESI.Blackbox((x1 + x2 - x1*x2)//(x1 + x2))
    vdhl = ESI.vanDerHoevenLecerf(R, 3, 2)
    (P, Q) = ESI.interpolate!(vdhl, f)
    @test P//Q == (x1 + x2 - x1*x2)//(x1 + x2)

    R, xi = PolynomialRing(QQ, ["x$i" for i in 1:10])
    cases = [
        sum(xi) // (prod(xi)),
        (xi[1] + xi[10])//(xi[1] - xi[10] - 2),
        (xi[1] + xi[3] + xi[5] + xi[7])^5//(xi[1] + 2),
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        vdhl = ESI.vanDerHoevenLecerf(R, n + 1, d + 1)
        P, Q = ESI.interpolate!(vdhl, f)
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

end

@testset "van-der-Hoeven-Lecerf, Z/Zp" begin
    R, (x1,x2) = PolynomialRing(GF(2^31-1), ["x1","x2"])
    cases = [
        (x1)//(2x2 + 1), R(2)//R(3), R(0)//R(1),
        (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        vdhl = ESI.vanDerHoevenLecerf(R, n, d)
        P, Q = ESI.interpolate!(vdhl, f)
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

    f = ESI.Blackbox((x1 + x2 - x1*x2)//(x1 + x2))
    vdhl = ESI.vanDerHoevenLecerf(R, 3, 2)
    (P, Q) = ESI.interpolate!(vdhl, f)
    @test P//Q == (x1 + x2 - x1*x2)//(x1 + x2)

    f = ESI.Blackbox((x1 - x2)^40//(2x1*x2))
    vdhl = ESI.vanDerHoevenLecerf(R, 40, 2)
    (P, Q) = ESI.interpolate!(vdhl, f)
    @test_broken P//Q == (x1 - x2)^40//(2x1*x2)

    R, xi = PolynomialRing(GF(fmpz(2)^300 + 157), ["x$i" for i in 1:10])
    cases = [
        (xi[1] - xi[2])^40//(2xi[1]*xi[2]),
        sum(xi) // (prod(xi)),
        (xi[1] + xi[10])//(xi[1] - xi[10] - 2),
        (xi[1] + xi[3] + xi[5] + xi[7])^5//(xi[1] + 2),
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        vdhl = ESI.vanDerHoevenLecerf(R, n + 1, d + 1)
        P, Q = ESI.interpolate!(vdhl, f)
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

end
