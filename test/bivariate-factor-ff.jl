ESI = ExactSparseInterpolations

@testset "factoring in Zp[x,y]" begin
    R, (x, y) = PolynomialRing(GF(2^31-1), ["x", "y"])
    
    cases = [
        (F = (y + 3)*(x + y)*((y + 1)*x + 2), Fi = [(y + 1)*x + 2, x + y]),

        (F = (x + 1)*(x + y), Fi = [x + 1, x + y],
        F2 = (x + 1)*(x + (3y)), Fi2 = [x + 1, x + (3y)]),
        
        (F = (y + 8)*(x + 1)*(x + y^4), Fi = [(x + 1), x + y^4],
        F2 = (x + 1)*(x + (88y)^4), Fi2 = [x + 1, x + (88y)^4]),
        
        (F = (((y + 1)^3 + 8)*x + 1)*(x^2 + y^4 + 1), Fi = [((y + 1)^3 + 8)*x + 1, x^2 + y^4 + 1],
        F2 = ((((2y) + 1)^3 + 8)*x + 1)*(x^2 + (2y)^4 + 1), Fi2 = [(((2y) + 1)^3 + 8)*x + 1, x^2 + (2y)^4 + 1]),

        (F = (x + y^20)*(x^20 + y + 1), Fi = [x + y^20, x^20 + y + 1]),

        (F = (x)*((y + 3)*x + 1)*((y + 4)*x + 2), Fi = [(x), ((y + 3)*x + 1), ((y + 4)*x + 2)]),
    ]

    for case in cases
        F, true_Fi = case.F, case.Fi
        if hasproperty(case, :F2)
            F2, true_Fi2 = case.F2, case.Fi2
        else
            F2, true_Fi2 = F, true_Fi
        end
        Fi, fi, Si = ESI.revealing_bivariate_factorization_ff(F)
        @test all(in(Fi), true_Fi) && all(in(true_Fi), Fi)
        @test Fi == ESI.bivariate_factorization_ff_with_known_univariate(F, fi, Si)
        Fi22 = ESI.bivariate_factorization_ff_with_known_univariate(F2, fi, Si)
        @test all(in(Fi22), true_Fi2) && all(in(true_Fi2), Fi22)
    end
end
