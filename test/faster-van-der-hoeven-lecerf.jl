ESI = ExactSparseInterpolations

@testset "fast vdH-Lecerf, Zp and Fq^n" begin
    for ground in [
        Nemo.GF(3*2^30+1), 
        Nemo.FqNmodFiniteField(fmpz(5), 12, :z)
    ]
        R, (x1,x2) = PolynomialRing(ground, ["x1","x2"])
        cases = [
            (x1)//(2x2 + 1), R(2)//R(3),
            (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3,
            (x1 + x2)^4//(x1 - x2)^6, (x1*x2 - x1 - x2)^3 // R(3)
        ]
        for case in cases
            f = ESI.Blackbox(case)
            @info "" f
            info = ESI.getboundsinfo(case)
            vdhl = ESI.FasterVanDerHoevenLecerf(R, info)
            P, Q = ESI.interpolate!(vdhl, f)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
        end

        R, xi = PolynomialRing(ground, ["x$i" for i in 1:10])
        cases = [
            sum(xi) // (prod(xi)),
            (xi[1] + xi[10])//(xi[1] - xi[10] - 2),
            (xi[1] + xi[3] + xi[5] + xi[7])^5//(xi[1] + 2),
        ]
        for case in cases
            f = ESI.Blackbox(case)
            info = ESI.getboundsinfo(case)
            vdhl = ESI.FasterVanDerHoevenLecerf(R, info)
            P, Q = ESI.interpolate!(vdhl, f)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
        end
    end
end
