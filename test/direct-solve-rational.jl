
ESI = ExactSparseInterpolations

@testset "Direct-Solve over rationals" begin
    for ground in [QQ, GF(2^31-1)]
        R, z = PolynomialRing(QQ, "z")
        cases = [
            (z^2)//(z + 1),
            R(1)//(z + 8),
            (z - 1)^10//(z + 1)^10,
            (z^3 - 3z^2 + 6z)//R(8),
            (z^3 - 3z^2 + 6z)//R(8),
            R(2)//R(3),
            ((z-1)(z-2)(z+3))//((z+6)^30)
        ]
        for case in cases
            f = ESI.Blackbox(case)
            n, d = degree(numerator(case)), degree(denominator(case))
            dsr = ESI.DirectSolveRational(R, n + 1, d + 1)
            P, Q = ESI.interpolate!(dsr, f)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
            dsr = ESI.DirectSolveRational(R, 2*n+1, 2*d+1)
            P, Q = ESI.interpolate!(dsr, f)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
        end
    end
end
