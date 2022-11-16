
ESI = ExactSparseInterpolations

@testset "Berlekamp-Massey over rationals" begin
    Ru, z = PolynomialRing(QQ, "z")
    Ring, Point, Poly = Nemo.FmpqPolyRing, Nemo.fmpq, Nemo.fmpq_poly
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

    bm = ESI.BerlekampMassey{Ring, Point, Poly}(Ru)

    for case in cases
        seq, steps, gen = case.seq, case.steps, case.gen
        Λ = nothing
        flag = nothing
        for (a, s) in zip(seq, steps)
            flag, Λ = ESI.next!(bm, Point(a))
            @test flag == s
        end
        if flag
            @test Λ == gen
        end
        empty!(bm)
    end
end

