ESI = ExactSparseInterpolations

@testset "fast gcd, Q and Z/Zp" begin
    for ground in [QQ, GF(2^31-1)]
        R, x = PolynomialRing(ground, "x")
        cases = [
            (R(0), R(1)), (R(4), R(0)),
            (R(4), R(2)), (x, R(3)), (x, x^2),
            (2x^3, 3x^4), (x + 1, 2(x + 1)),
            (x^2 + 1, x^2 + 2), (x^2 + 1, x^2 + 3),
            ((x - 1)*(x + 4)*(x - 8), (x - 8)*(x + 11)),
            ((x^4 + 4)*(x^2 + 3)*(x + 2), (x^2 + 2)*(x^4 + 4)),
            ((x - 2)^10*(x^4 + 4)^5*(x + 11), (x - 2)^7*(x^2 + 2)^11*(x + 11)),
            ((x - 1)^1000*(x - 2)^200, (x - 1)^500*(x - 2)^700)
        ]
        for (P, Q) in cases
            @test Nemo.gcd(P, Q) == ESI.fastgcd(P, Q)
        end
    
        # # random tests
        for i in 1:100
            getrandpoly = () -> rand(R, 1:50, 1:5)
            a,b,c = getrandpoly(), getrandpoly(), getrandpoly()
            P = (1 + a*rand(0:1))*(1 + b*rand(0:1))*(1 + c*rand(0:1))
            Q = (1 + a*rand(0:1))*(1 + b*rand(0:1))*(1 + c*rand(0:1))
            @test Nemo.gcd(P, Q) == ESI.fastgcd(P, Q)
        end
    end
end
