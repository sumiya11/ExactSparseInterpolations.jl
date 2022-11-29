ESI = ExactSparseInterpolations

@testset "v-d-h-l-GCD, Q" begin
    R, (x1,x2) = PolynomialRing(QQ, ["x1","x2"])
    cases = [
        (R(1), R(2)), (x1, R(3)), (x1*x2, x1^2),
        ((x1 + x2)^3, (x1 + x2)^4), ((x1 + x2), 2(x1 + x2)),
        ((x1*x2 - x1)*(x1 + x2)^3*(x2 - x1)^4, (x1*x2 - x1)*(x2 - x1)^3)
    ]
    for (P, Q) in cases
        f1, f2 = ESI.Blackbox(P), ESI.Blackbox(Q)
        n, d = maxdeg(P), maxdeg(Q)
        vdhlg = ESI.vanDerHoevenLecerfGCD(R, n, d)
        g = ESI.gcd!(vdhlg, f1, f2)
        @test isone(leading_coefficient(g))
        @test Nemo.gcd(P, Q) == g
    end

    # random tests
    R, (x1,x2,x3,x4) = PolynomialRing(QQ, ["x1","x2","x3","x4"])
    for i in 1:50
        getrandpoly = () -> rand(R, 1:2, 1:2, 1:5)
        P = getrandpoly()*getrandpoly()*getrandpoly()
        Q = getrandpoly()*getrandpoly()*getrandpoly()
        f1, f2 = ESI.Blackbox(P), ESI.Blackbox(Q)
        n, d = maxdeg(P), maxdeg(Q)
        vdhlg = ESI.vanDerHoevenLecerfGCD(R, n, d)
        g = ESI.gcd!(vdhlg, f1, f2)
        @test isone(leading_coefficient(g))
        @test Nemo.gcd(P, Q) == g
    end
end
