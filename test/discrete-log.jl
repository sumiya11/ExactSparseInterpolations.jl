ESI = ExactSparseInterpolations

# good prime ? 2 * 3 * 7 * 47 * 769 * 1193 * 4583 * 8167 * 17417

@testset "baby-step-giant-step, Pohlig Hellman" begin
    F = Nemo.GF(17)
    PF = ESI.PrecomputedField(F)
    buf = ESI.DiscreteLogBuffers(PF)
    ord = 16
    # for each generator of F
    for a in [3, 5, 6, 7, 10, 11, 12, 14]
        a = F(a)
        for d in 0:15
            y = a^d
            x = ESI.direct_discrete_log(a, y, ord, buf)
            @test x == d
            x = ESI.babystep_giantstep_discrete_log(a, y, ord, buf)
            @test x == d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, buf)
            @test x == d
        end
    end
    F = Nemo.GF(101)
    ord = 100
    PF = ESI.PrecomputedField(F)
    buf = ESI.DiscreteLogBuffers(PF)
    # for some generators of F
    for a in [2, 3, 7, 8, 11, 12, 15, 18, 26, 27, 28, 29, 34, 35, 38, 40, 42, 46, 48, 50, 98, 99]
        a = F(a)
        for d in 0:99
            y = a^d
            x = ESI.direct_discrete_log(a, y, ord, buf)
            @test x == d
            x = ESI.babystep_giantstep_discrete_log(a, y, ord, buf)
            @test x == d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, buf)
            @test x == d
        end
    end
    F = Nemo.GF(2^16+1)
    ord = 2^16
    PF = ESI.PrecomputedField(F)
    buf = ESI.DiscreteLogBuffers(PF)
    # for some generators of F
    for a in [3, 5, 6, 7, 10, 11, 12, 14, 20, 22, 23, 24]
        a = F(a)
        for d in 0:2^10:2^16-1
            y = a^d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, buf)
            @test x == d
            x = ESI.discrete_log(a, y, ord, buf)
            @test x == d
        end
    end
    F = Nemo.GF(3*2^30+1)
    ord = 3*2^30
    PF = ESI.PrecomputedField(F)
    buf = ESI.DiscreteLogBuffers(PF)
    # for some generators of F
    for a in [274810138, 2082094268, 1727712020, 1894941899]
        a = F(a)
        for d in rand(0:3*2^30-1, 100)
            y = a^d
            x = ESI.pohlig_hellman_discrete_log(a, y, ord, buf)
            @test x == d
            x = ESI.discrete_log(a, y, ord, buf)
            @test x == d
        end
    end
end

@testset "discrete log, base isn't a generator" begin
    F = Nemo.GF(17)
    ord = 16
    PF = ESI.PrecomputedField(F)
    buf = ESI.DiscreteLogBuffers(PF)
    # for each non-generator of F
    for a in [2, 4, 8, 9, 13, 15]
        a = F(a)
        a_ord = findfirst(map(i -> isone(a^i), 1:ord))
        for d in 0:a_ord-1
            y = a^d
            x = ESI.direct_discrete_log(a, y, a_ord, buf)
            @test x == d
            x = ESI.babystep_giantstep_discrete_log(a, y, a_ord, buf)
            @test x == d
            x = ESI.pohlig_hellman_discrete_log(a, y, a_ord, buf)
            @test x == d
            x = ESI.discrete_log(a, y, a_ord, buf)
            @test x == d
        end
    end

    F = Nemo.GF(101)
    PF = ESI.PrecomputedField(F)
    buf = ESI.DiscreteLogBuffers(PF)
    g = F(83) # g is a generator
    for (a, a_ord) in [(g^5, 20), (g^4, 25), (g^2, 50), (g^10, 10)]
        a = F(a)
        for d in 0:a_ord-1
            y = a^d
            x = ESI.direct_discrete_log(a, y, a_ord, buf)
            @test x == d
            x = ESI.babystep_giantstep_discrete_log(a, y, a_ord, buf)
            @test x == d
            x = ESI.pohlig_hellman_discrete_log(a, y, a_ord, buf)
            @test x == d
            x = ESI.discrete_log(a, y, a_ord, buf)
            @test x == d
        end
    end
end
