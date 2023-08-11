ESI = ExactSparseInterpolations

@testset "univariate Ben-or-Tiwari, Zp" begin
    for ground in [
        Nemo.GF(3 * 2^30 + 1),
    ]
        R, (x1,) = PolynomialRing(ground, ["x1"])
        cases = [
            R(1), x1, 2 * x1, x1^3 + x1,
            x1^10 + 1, (x1 + 2)^6,
            x1^10000 + 2x1^2000 - 100x1^2 + 1,
            (x1 + 1)^100, (x1 + 1)^500
        ]
        for case in cases
            T = length(case)
            f = ESI.Blackbox(case)

            bot = ESI.PrimesBenOrTiwari(R, T)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1)
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case
        end
    end
end

@testset "primes Ben-or-Tiwari, Zp" begin
    for ground in [
        Nemo.GF(3 * 2^30 + 1),
    ]
        R, (x1, x2, x3) = PolynomialRing(ground, ["x1", "x2", "x3"])
        cases = [
            R(1), x1, 2 * x1, x1^3 + x1,
            (x1 + x2)^10 + 1, (x1 + 2)^6 - x2^10,
            x1^10000 + 2x2^2000 - 100x3^2 + 1,
            (x1 + 1)^100,
            (x1 + x2 + x3 + 1)^10,
            (x1 + x2^2 + x3^5 + 1)^10,
            x1^(2^9) + 2x2^(2^7) + 3x3^(2^10),
            12321(x1 * x2)^10, (x1 * x2 * x3)^40,
            x1 * x2, x2 * x3, x1 * x3
        ]
        for case in cases
            T = length(case)
            f = ESI.Blackbox(case)
            pd = map(x -> degree(case, x), gens(R))

            bot = ESI.PrimesBenOrTiwari(R, T)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1)
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case

            # Test for the case when pd is an upper bound
            as = rand(0:5, length(pd))
            pd = pd .+ as
            bot = ESI.PrimesBenOrTiwari(R, T)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1)
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case
        end
    end
end

@testset "primes Ben-or-Tiwari, upper bound on T" begin
    for ground in [
        Nemo.GF(3 * 2^30 + 1),
    ]
        R, (x1,) = PolynomialRing(ground, ["x1"])
        cases = [
            R(1), x1, 2 * x1, x1^3 + x1,
            x1^10 + 1, (x1 + 2)^6,
            x1^10000 + 2x1^2000 - 100x1^2 + 1,
            (x1 + 1)^100, (x1 + 1)^500
        ]
        for case in cases
            T = length(case)
            f = ESI.Blackbox(case)

            a = rand(1:100)
            T = a * T
            bot = ESI.PrimesBenOrTiwari(R, T)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1)
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case
        end

        R, (x1, x2, x3) = PolynomialRing(ground, ["x1", "x2", "x3"])
        cases = [
            R(1), x1, 2 * x1, x1^3 + x1,
            (x1 + x2)^10 + 1, (x1 + 2)^6 - x2^10,
            x1^10000 + 2x2^2000 - 100x3^2 + 1,
            (x1 + 1)^100
        ]
        for case in cases
            T = length(case)
            f = ESI.Blackbox(case)
            pd = map(x -> degree(case, x), gens(R))

            a = rand(1:100)
            T = a * T
            bot = ESI.FasterMultivariateBenOrTiwari(R, T, pd)
            ω = ESI.startingpoint(bot)
            ωs = map(i -> ω .^ i, 0:2T-1)
            ys = map(f, ωs)
            @test ESI.interpolate!(bot, ωs, ys) == case
        end
    end
end
