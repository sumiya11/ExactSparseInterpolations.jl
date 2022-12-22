ESI = ExactSparseInterpolations

@testset "fast Cauchy, Q and Z/Zp" begin
    for ground in [QQ, GF(2^31-1)]
        R, z = PolynomialRing(ground, "z")
        cases = [
            (R(0))//(R(1)),
            (z^2)//(z + 1),
            R(1)//(z + 8),
            (z - 1)^5//(z + 1)^5,
            (z^3 - 3z^2 + 6z)//R(8),
            (z^3 - 3z^2 + 6z)//R(8),
            R(2)//R(3),
            ((z-1)(z-2)(z+3))//((z+6)^5),
            (3z - 1)//R(4),
            (z^8 + 1)//(z^4 + 1),
            (z^16 + 1)//(z^8 + 1),
            (z + 5)^10//(z - 1)^7,
        ]
        for case in cases
            f = ESI.Blackbox(case)
            n, d = max(0, degree(numerator(case))), degree(denominator(case))
            c = ESI.FasterCauchy(R, n, d)
            xs = ESI.distinct_points(ground, n + d + 2)
            ys = map(f, xs)
            P, Q = ESI.interpolate!(c, xs, ys)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case

            # Test for the case when n and d are upper bounds
            a, b = rand(1:10), rand(1:10)
            n, d = a*max(0, degree(numerator(case))), b*max(0, degree(denominator(case)))
            c = ESI.FasterCauchy(R, n, d)
            xs = ESI.distinct_points(ground, n + d + 2)
            ys = map(f, xs)
            P, Q = ESI.interpolate!(c, xs, ys)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
        end
    end
end
