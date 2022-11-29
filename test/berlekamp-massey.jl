
ESI = ExactSparseInterpolations

@testset "Berlekamp-Massey, Q" begin
    K = QQ
    Ru, z = PolynomialRing(K, "z")
    cases = [
        (
            seq=[1, 2, 4, 8, 16], 
            steps=[false, false, true, true, true], 
            gen=-2z + 1
        ),
        (
            seq=[0, 1, 1, 2, 3, 5, 8, 13], 
            steps=[true, false, false, true, true, true, true, true], 
            gen=-z^2 - z + 1
        ),
        (
            seq=[0, 1, 4, 9, 16, 25, 36, 49], 
            steps=[true, false, false, false, false, false, true, true], 
            gen=-z^3 + 3*z^2 - 3*z + 1
        ),
        (
            seq=[1, 1, 1, 1, 1], 
            steps=[false, true, true, true, true], 
            gen=-z + 1
        ),
        (
            seq=[1, 1, 2, 5, 14, 42, 132, 429], 
            steps=[false, true, false, false, false, false, false, false], 
            gen=-z + 1
        ),
    ]

    bm = ESI.BerlekampMassey(Ru)

    for case in cases
        seq, steps, gen = case.seq, case.steps, case.gen
        Λ = nothing
        flag = nothing
        for (a, s) in zip(seq, steps)
            flag, Λ = ESI.next!(bm, K(a))
            @test flag == s
        end
        if flag
            @test Λ == gen
        end
        empty!(bm)
    end
end

@testset "Berlekamp-Massey, Z/Zp" begin
    K = GF(5)
    Ru, z = PolynomialRing(K, "z")
    cases = [
        (
            seq=[1, 2, 4, 8, 16], 
            steps=[false, false, true, true, true], 
            gen=-2z + 1
        ),
        (
            seq=[0, 1, 1, 2, 3, 5, 8, 13], 
            steps=[true, false, false, true, true, true, true, true], 
            gen=-z^2 - z + 1
        ),
        (
            seq=[1, 1, 1, 1, 1], 
            steps=[false, true, true, true, true], 
            gen=4z + 1
        ),
        (
            seq=[2, 1, 3, 4, 2],
            steps=[false, false, true, true, true], 
            gen=2z + 1
        ),
    ]

    for case in cases
        bm = ESI.BerlekampMassey(Ru)
        seq, steps, gen = case.seq, case.steps, case.gen
        Λ = nothing
        flag = nothing
        for (a, s) in zip(seq, steps)
            flag, Λ = ESI.next!(bm, K(a))
            @test flag == s
        end
        if flag
            @test Λ == gen
        end
    end

end
