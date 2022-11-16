
ESI = ExactSparseInterpolations

@testset "Newton over rationals" begin
    R, x1 = PolynomialRing(QQ, "x1")
    cases = [
        R(1), x1, 2*x1, x1^3 + x1, 
        x1^10 + 1, (x1 + 2)^6
    ]
    for case in cases
        f = ESI.Blackbox(case)
        N = ESI.Newton(R, d=degree(case)) 
        @test ESI.interpolate!(N, f) == case
        N = ESI.Newton(R) 
        @test ESI.interpolate!(N, f) == case
    end
    for case in cases
        f = ESI.Blackbox(case)
        N = ESI.Newton(R)
        for d in 0:degree(case)-1
            x = ESI.random_point(QQ)
            y = f(x)
            success, poly = ESI.next!(N, x, y)
        end
        x = ESI.random_point(QQ)
        y = f(x)
        @test ESI.next!(N, x, y) == (false, case)
        x = ESI.random_point(QQ)
        y = f(x)
        @test ESI.next!(N, x, y) == (true, case)
    end
end
