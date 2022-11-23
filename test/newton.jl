
ESI = ExactSparseInterpolations

@testset "Newton, Q and Z/Zp" begin
    for ground in [QQ, GF(2^31-1)]
        R, x1 = PolynomialRing(ground, "x1")
        cases = [
            R(1), x1, 2*x1, x1^3 + x1, 
            x1^10 + 1, (x1 + 2)^6,
            (x1 - 8)^10*(4x1 + 1)^5*(x1^3 + 2x1 + 1)^10
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
                x = ESI.random_point(ground)
                y = f(x)
                success, poly = ESI.next!(N, x, y)
            end
            x = ESI.random_point(ground)
            y = f(x)
            @test ESI.next!(N, x, y) == (false, case)
            x = ESI.random_point(ground)
            y = f(x)
            @test ESI.next!(N, x, y) == (true, case)
        end
    end
end
