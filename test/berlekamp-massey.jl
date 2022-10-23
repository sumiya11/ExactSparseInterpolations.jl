
EI = ExactInterpolations

@testset "Berlekamp Massey over rationals" begin
    Ru, z = PolynomialRing(QQ, "z")
    I, R, Y = Nemo.fmpq, Nemo.FmpqPolyRing, Nemo.fmpq_poly
    cases = [
        (
            seq=[1, 2, 4, 8, 16], 
            steps=[false, false, true, true, true], 
            gen=-2z + 1
        ),
        (
            seq=[0, 1, 1, 2, 3, 5, 8, 13], 
            steps=[true, false, false, false, true, true, true, true], 
            gen=-z^2 - z + 1
        ),
        (
            seq=[0, 1, 4, 9, 16, 25, 36, 49], 
            steps=[true, false, false, false, false, false, true, true], 
            gen=-z^3 + 3*z^2 - 3*z + 1
        ),
        (
            seq=[1, 1, 1, 1, 1], 
            steps=[false, false, true, true, true], 
            gen=-z + 1
        ),
        (
            seq=[1, 1, 2, 5, 14, 42, 132, 429], 
            steps=[false, false, false, false, false, false, false, false], 
            gen=-z + 1
        ),
    ]

    bm = EI.BerlekampMassey{I, R, Y}(Ru)

    for case in cases
        seq, steps, gen = case.seq, case.steps, case.gen
        Λ = nothing
        flag = nothing
        for (a, s) in zip(seq, steps)
            flag, Λ = EI.next!(bm, I(a))
            @test flag == s
        end
        if flag
            @test Λ == gen
        end
        empty!(bm)
    end
end
