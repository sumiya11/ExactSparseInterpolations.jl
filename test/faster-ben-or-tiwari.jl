ESI = ExactSparseInterpolations

@testset "fast univ. Ben-or-Tiwari, Zp and Fq^n" begin
    for ground in [
        Nemo.GF(3*2^30+1), 
        Nemo.FqNmodFiniteField(fmpz(5), 6, :z)
    ]
        R, (x1,) = PolynomialRing(ground, ["x1"])
        cases = [
            R(1), x1, 2*x1, x1^3 + x1, 
            x1^10 + 1, (x1 + 2)^6,
            x1^10000 + 2x1^2000 - 100x1^2 + 1,
            (x1 + 1)^100, (x1 + 1)^1000
        ]
        for case in cases
            T = length(case)
            f = ESI.Blackbox(case)
            
            bot = ESI.FasterBenOrTiwari(R, T)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1) 
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case

            bot = ESI.FasterBenOrTiwari(R, T)
            @test ESI.interpolate!(bot, f) == case
        end
    end
end

@testset "fast mult. Ben-or-Tiwari, Zp and Fq^n" begin
    for ground in [
        Nemo.GF(3*2^30+1), 
        Nemo.FqNmodFiniteField(fmpz(5), 12, :z)
    ]
        R, (x1,x2,x3) = PolynomialRing(ground, ["x1","x2","x3"])
        cases = [
            R(1), x1, 2*x1, x1^3 + x1, 
            (x1 + x2)^10 + 1, (x1 + 2)^6 - x2^10,
            x1^10000 + 2x2^2000 - 100x3^2 + 1,
            (x1 + 1)^100,
            (x1 + x2 + x3 + 1)^10,
            (x1 + x2^2 + x3^5 + 1)^10,
            x1^(2^9) + 2x2^(2^7) + 3x3^(2^10),
            12321(x1*x2)^10, (x1*x2*x3)^40
        ]
        for case in cases
            T = length(case)
            f = ESI.Blackbox(case)
            pd = map(x -> degree(case, x), gens(R))

            bot = ESI.FasterMultivariateBenOrTiwari(R, T, pd)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1) 
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case

            bot = ESI.FasterMultivariateBenOrTiwari(R, T, pd)
            @test ESI.interpolate!(bot, f) == case
        end
    end
end
