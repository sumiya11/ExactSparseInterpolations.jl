
ESI = ExactSparseInterpolations

# degrees over 30 over rationals are quite hard

@testset "Cauchy, Q and Z/Zp" begin
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
            (3z - 1)//R(4)
        ]
        for case in cases
            f = ESI.Blackbox(case)
            n, d = degree(numerator(case)), degree(denominator(case))
            c = ESI.Cauchy(R, n + 1, d + 1)
            P, Q = ESI.interpolate!(c, f)
            @test isone(leading_coefficient(Q))
            @test P//Q == case
            c = ESI.Cauchy(R, 2*(n+1), 2*(d+1))
            P, Q = ESI.interpolate!(c, f)
            @test isone(leading_coefficient(Q))
            @test P//Q == case
        end
    end
end
