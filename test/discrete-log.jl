ESI = ExactSparseInterpolations

# good prime 2 * 3 * 7 * 47 * 769 * 1193 * 4583 * 8167 * 17417

@testset "Discrete logarithms" begin
    F = Nemo.GF(17)
    ord = 16
    factors = collect(Primes.factor(Dict, ord))
    # for each generator of F
    for a in [3, 5, 6, 7, 10, 11, 12, 14]
        a = F(a)
        for d in 0:15
            y = a^d
            x = ESI.direct_discrete_log(a, y)
            @test x == d
            x = ESI.babystep_giantstep_discrete_log(a, y, ord)
            @test x == d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, factors)
            @test x == d
        end
    end
    F = Nemo.GF(101)
    ord = 100
    factors = collect(Primes.factor(Dict, ord))
    # for some generators of F
    for a in [2, 3, 7, 8, 11, 12, 15, 18, 26, 27, 28, 29, 34, 35, 38, 40, 42, 46, 48, 50, 98, 99]
        a = F(a)
        for d in 0:99
            y = a^d
            x = ESI.direct_discrete_log(a, y)
            @test x == d
            x = ESI.babystep_giantstep_discrete_log(a, y, ord)
            @test x == d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, factors)
            @test x == d
        end
    end
    F = Nemo.GF(2^16+1)
    ord = 2^16
    factors = collect(Primes.factor(Dict, ord))
    # for some generators of F
    for a in [3, 5, 6, 7, 10, 11, 12, 14, 20, 22, 23, 24]
        a = F(a)
        for d in 0:2^10:2^16-1
            y = a^d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, factors)
            @test x == d
        end
    end
end
