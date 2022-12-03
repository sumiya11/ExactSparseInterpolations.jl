
maxdeg(f) = maximum(map(x -> degree(f, x), gens(parent(f)))) 

ESI = ExactSparseInterpolations

@testset "Ben-or-Tiwari, Q" begin
    R, (x1,) = PolynomialRing(QQ, ["x1"])
    cases = [
        R(1), x1, 2*x1, x1^3 + x1, 
        x1^10 + 1, (x1 + 2)^6
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end

    R, (x1,x2) = PolynomialRing(QQ, ["x1","x2"])
    cases = [
        R(1), x1 + x2, x1^2 + x2^2, x2,
        (x1 + 2x2)^10, x1^100, x2^100
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end

    R, x = PolynomialRing(QQ, ["x$i" for i in 1:10])
    cases = [
        R(1), sum(x), prod(x),
        (x[1] + 2*x[3] + 3*x[5] + 4*x[7] + 5*x[9])^2,
        x[1]^20 + x[10]^30
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end
end

@testset "Ben-or-Tiwari, Z/Zp" begin
    R, (x1,) = PolynomialRing(GF(2^31-1), ["x1"])
    cases = [
        R(1), x1, 2*x1, x1^3 + x1, 
        x1^10 + 1, (x1 + 2)^6
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end

    R, (x1,x2) = PolynomialRing(GF(2^31-1), ["x1","x2"])
    cases = [
        R(1), x1 + x2, x1^2 + x2^2, x2,
        (x1 + 2x2)^5
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end

    R, x = PolynomialRing(GF(2^31-1), ["x$i" for i in 1:10])
    cases = [
        R(1), sum(x),
        x[1]^25 + x[2]^2
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end

    R, (x,y,z) = PolynomialRing(GF(Nemo.fmpz(2)^300 + 157), ["x","y","z"])
    cases = [
        x^10 + y^10 + z^10,
        x^100 + (y - 1)^50 + (z - 1)^30
    ]
    for case in cases
        f = ESI.Blackbox(case)
        Z = ESI.BenOrTiwari(R) 
        @test ESI.interpolate!(Z, f) == case
    end
end
