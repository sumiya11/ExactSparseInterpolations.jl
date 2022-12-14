ESI = ExactSparseInterpolations

@testset "simultaneous adaptive v-d-h-l" begin
    R1, (x1,x2,x3) = PolynomialRing(GF(2^31-1), ["x1","x2","x3"])
    cases = [
        [3x1//R1(1)],
        [R1(1)//R1(x2^3), R1(1)//(x2), (x1*x3)//(x2 + 8)],
        [R1(1)//R1(1), R1(1)//(x2^4), (x1 - x2)^2 // (x2 - x3)^2],
        [R1(4)//R1(5), x2^4//x1^3],
        [(x1 + x2 + x3)^10//(x1 - x2 - 1)^15, R1(1)//x1^16],
        [x1^100 // (x2^200 + 1), x1^200 // (x3^100 + 2)],
        [(x1^1001 + x2^2001 + x3^3001) // (x1 + x2 + x3)]
    ]
    for case in cases
        fs = [ESI.Blackbox(c) for c in case]
        count = length(case)
        N = maximum(map(total_degree ∘ numerator, case))
        D = maximum(map(total_degree ∘ denominator, case))
        Is = ESI.SimultaneousAdaptiveVanDerHoevenLecerf(R1, count, N, D)
        while !ESI.allready(Is)
            xs = ESI.nextpoints!(Is)
            ys = [f.(xs) for f in fs]
            ESI.nextevaluations!(Is, ys)
        end
        for i in 1:count
            ans = numerator(case[i]), denominator(case[i])
            @test ESI.getresult(Is, i) == (true, ans...)
        end
    end

    R2, ys = PolynomialRing(GF(2^31-1), ["y$i" for i in 1:31])
    cases = [
        [sum(ys) // prod(ys), prod(ys) // sum(ys)],
        [R2(1) // sum(ys[i]^i for i in 1:10), (ys[1] + ys[2] + ys[3])^10 // (ys[31]^13 + 111)]
    ]
    for case in cases
        fs = [ESI.Blackbox(c) for c in case]
        count = length(case)
        N = maximum(map(total_degree ∘ numerator, case))
        D = maximum(map(total_degree ∘ denominator, case))
        Is = ESI.SimultaneousAdaptiveVanDerHoevenLecerf(R2, count, N, D)
        while !ESI.allready(Is)
            xs = ESI.nextpoints!(Is)
            ys = [f.(xs) for f in fs]
            ESI.nextevaluations!(Is, ys)
        end
        for i in 1:count
            ans = numerator(case[i]), denominator(case[i])
            @test ESI.getresult(Is, i) == (true, ans...)
        end
    end
end
