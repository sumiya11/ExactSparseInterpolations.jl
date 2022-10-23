
EI = ExactInterpolations

@testset "Newton over rationals" begin
    R, x1 = PolynomialRing(QQ, "x1")
    cases = [
        R(1), x1, 2*x1, x1^3 + x1, 
        x1^10 + 1, (x1 + 2)^6
    ]
    for case in cases
        f = EI.Blackbox(case)
        N = EI.Newton(R, d=degree(case)) 
        @test interpolate(N, f) == case
    end
end
