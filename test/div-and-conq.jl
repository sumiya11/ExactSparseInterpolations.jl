ESI = ExactSparseInterpolations

@testset "fast univariate interpolate, Z/Zp, Fq^n" begin
    for ground in [
            Nemo.GF(2^31-1), 
            Nemo.FqNmodFiniteField(fmpz(5), 5, :z)
        ]
        R, x = PolynomialRing(ground, "x")
        cases = [
            R(0), R(4), x, x^2, 2x^3, 3x^4, x + 2,
            x^2 + 1, (x - 1)*(x + 4)*(x - 8),
            (x^4 + 4)*(x^2 + 3)*(x + 2),
            (x - 2)^10*(x^4 + 4)^5*(x + 11),
            (x - 1)^10*(x - 2)^20, (x - 1)^50*(x - 2)^70,
            x^20 + 1,
        ]
        for f in cases
            xs = ESI.distinct_points(ground, max(degree(f) + 1, 1))
            ys = map(x -> evaluate(f, x), xs)
            @test f == Nemo.interpolate(R, xs, ys) == ESI.fastpolyinterpolate(R, xs, ys)
        end
    
        # random tests
        for i in 1:100
            getrandpoly = () -> rand(R, 1:50, 1:5)
            a,b,c = getrandpoly(), getrandpoly(), getrandpoly()
            f = (1 + a*rand(0:1))*(1 + b*rand(0:1))*(1 + c*rand(0:1))
            xs = ESI.distinct_points(ground, max(degree(f) + 1, 1))
            ys = map(x -> evaluate(f, x), xs)
            @test f == Nemo.interpolate(R, xs, ys) == ESI.fastpolyinterpolate(R, xs, ys)
        end
    end
end

@testset "fast vandermonde solve, Z/Zp, Fq^n" begin
    for ground in [
            Nemo.GF(2^31-1), 
            Nemo.FqNmodFiniteField(fmpz(5), 5, :z)
        ]
        R, (x,) = PolynomialRing(ground, ["x"])
        Runiv, xuniv = PolynomialRing(ground, "x")
        cases = [
            R(4), x, x^2, 2x^3, 3x^4, x + 2,
            x^2 + 1, (x - 1)*(x + 4)*(x - 8),
            (x^4 + 4)*(x^2 + 3)*(x + 2),
            (x - 2)^10*(x^4 + 4)^5*(x + 11),
            (x - 1)^10*(x - 2)^20, (x - 1)^50*(x - 2)^70,
            x^20 + 1, (x - 1)^100
        ]
        for f in cases
            x = [ESI.randomgenerator(ground)]
            vs = map(m -> evaluate(m, x), collect(monomials(f)))
            ys = map(i -> evaluate(f, x .^ i), 0:length(vs)-1)
            @test collect(coefficients(f)) == ESI.solve_transposed_vandermonde(Runiv, vs, ys)
        end
    end
end
