
maxdeg(f) = maximum(map(x -> degree(f, x), gens(parent(f)))) 

EI = ExactSparseInterpolations

@testset "Zippel over rationals" begin
    R, (x1,) = PolynomialRing(QQ, ["x1"])
    cases = [
        R(1), x1, 2*x1, x1^3 + x1, 
        x1^10 + 1, (x1 + 2)^6
    ]
    for case in cases
        f = EI.Blackbox(case)
        Z = EI.Zippel(R, maxdeg(case)+1) 
        @test interpolate(Z, f) == case
    end

    R, (x1,x2) = PolynomialRing(QQ, ["x1","x2"])
    cases = [
        R(1), x1 + x2, x1^2 + x2^2,
        (x1 + 2x2)^10, x1^100, x2^100
    ]
    for case in cases
        f = EI.Blackbox(case)
        Z = EI.Zippel(R, maxdeg(case)+1) 
        @test interpolate(Z, f) == case
    end

    R, x = PolynomialRing(QQ, ["x$i" for i in 1:10])
    cases = [
        R(1), sum(x), prod(x),
        (x[1] + 2*x[3] + 3*x[5] + 4*x[7] + 5*x[9])^2,
        x[1]^20 + x[10]^30
    ]
    for case in cases
        f = EI.Blackbox(case)
        Z = EI.Zippel(R, maxdeg(case)+1) 
        @test interpolate(Z, f) == case
    end
end
