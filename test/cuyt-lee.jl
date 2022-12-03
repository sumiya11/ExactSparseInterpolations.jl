ESI = ExactSparseInterpolations

maxdeg(f) = maximum(map(t -> total_degree(t), terms(f)), init=0) 

@testset "Cuyt-Lee, Q" begin
    R, (x1,x2) = PolynomialRing(QQ, ["x1","x2"])
    cases = [
        (x1)//(2x2 + 1), R(2)//R(3), R(0)//R(1),
        (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        cl = ESI.CuytLee(R, n + 1, d + 1)
        P, Q = ESI.interpolate!(cl, f)
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

    R, xi = PolynomialRing(QQ, ["x$i" for i in 1:10])
    cases = [
        sum(xi) // (prod(xi)),
        (xi[1] + xi[10])//(xi[1] - xi[10] - 2)
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        cl = ESI.CuytLee(R, n + 1, d + 1)
        P, Q = ESI.interpolate!(cl, f)
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

    f = ESI.Blackbox((x1 + x2 - x1*x2)//(x1 + x2))
    cl = ESI.CuytLee(R, 3, 2)
    @test_broken (P, Q) = ESI.interpolate!(cl, f)
end
