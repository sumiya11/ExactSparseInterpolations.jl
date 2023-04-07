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

        # for t == T
        # so that the number of points equals the number of terms
        cases = [
            R(4), x, x^2, 2x^3, 3x^4, x + 2,
            x^2 + 1, (x - 1)*(x + 4)*(x - 8),
            (x^4 + 4)*(x^2 + 3)*(x + 2),
            (x - 2)^10*(x^4 + 4)^5*(x + 11),
            (x - 1)^10*(x - 2)^20, (x - 1)^50*(x - 2)^70,
            x^20 + 1, (x - 1)^100
        ]
        for f in cases
            ω = [ESI.randomgenerator(ground)]
            vs = map(m -> evaluate(m, ω), collect(monomials(f)))
            ys = map(i -> evaluate(f, ω .^ i), 0:length(vs)-1)
            @test collect(coefficients(f)) == ESI.solve_transposed_vandermonde(Runiv, vs, ys)
        end

    end

end

@testset "fast toeplitz solve, Z/Zp" begin
    K = GF(2^31-1)
    R, t = K["t"]

    a = map(K, [1, 1, 1, 2, 3, 4, 5])
    b = map(K, [5, 7, 10, 14])
    @test t^3 + t^2 + t + 1 == ESI.solve_toeplitz(R, a, b)
    a = map(K, [1, 1, 2])
    b = map(K, [5, 7])
    @test 3t + 2 == ESI.solve_toeplitz(R, a, b)
    a = map(K, [1])
    b = map(K, [5])
    @test 5 == ESI.solve_toeplitz(R, a, b)

    for K in [
        Nemo.GF(2^31-1), 
        Nemo.FqNmodFiniteField(fmpz(5), 5, :z)
    ]
        R, t = PolynomialRing(K, "t")
        for n in (2, 3, 5, 10, 100, 300)
            b = map(_ -> rand(K), 1:n)
            a = map(_ -> rand(K), 1:2n-1)
            S = MatrixSpace(K, n, n)
            s = MatrixSpace(K, n, 1)
            A = zero(S)
            for i in 1:2n-1
                if i <= n
                    for j in n - i + 1:-1:1
                        A[j, j + i - 1] = a[n - i + 1]
                    end
                else
                    for j in (i - n + 1):n
                        A[j, j - (i - n)] = a[i]
                    end
                end
            end
            iszero(det(A)) && continue
            solution_1 = collect(coefficients(ESI.solve_toeplitz(R, a, b)))
            solution_true = collect(Nemo.solve(A, s(b)))[:, 1]
            @test solution_1 == solution_true
        end
    end
end

@testset "fast multivariate evaluate, Z/Zp" begin
    K = Nemo.GF(2^31-1)
    R, (x,y) = PolynomialRing(K, ["x","y"])

    T = 3
    f = x + y^2 + 3
    ω = [K(2), K(3)]
    @test [K(5), K(14), K(88)] == ESI._fast_multivariate_evaluate(R, f, ω, T)

    for K in [
            GF(2^31-1), 
            GF(4611686018427388039)
            # FqNmodFiniteField(fmpz(7), 5, :z)
        ]
        R, (x,y,z) = PolynomialRing(K, ["x","y","z"])
        cases = [
            (poly=x*y*z + y + z, ω=[K(2),K(3),K(5)]),
            (poly=x + y + z, ω=[K(2),K(3),K(5)])
        ]
        for case in cases
            poly, ω = case.poly, case.ω
            t = length(poly)
            evals_1 = ESI._fast_multivariate_evaluate(R, poly, ω, t)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:t-1)
            @test evals_1 == evals_true

            k = rand(t+1:10t)
            evals_1 = ESI._fast_multivariate_evaluate(R, poly, ω, k)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:k-1)
            @test evals_1 == evals_true
        end
        ω = [K(5), K(3), K(2)]
        for i in 1:100
            poly = rand(R, 1:100, 1:10)
            t = length(poly)
            evals_1 = ESI._fast_multivariate_evaluate(R, poly, ω, t)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:t-1)
            @test evals_1 == evals_true
        end
        R, xi = PolynomialRing(K, ["x$i" for i in 1:7])
        ω = [K(5), K(3), K(2), K(17), K(11), K(7), K(13)]
        for i in 1:10
            poly = rand(R, 1:1000, 0:10)
            t = length(poly)
            evals_1 = ESI._fast_multivariate_evaluate(R, poly, ω, t)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:t-1)
            @test evals_1 == evals_true

            k = rand(t+1:10t)
            evals_1 = ESI._fast_multivariate_evaluate(R, poly, ω, k)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:k-1)
            @test evals_1 == evals_true
        end
    end

    K = Nemo.GF(2^31-1)
    R, (x,y) = PolynomialRing(K, ["x","y"])
    f = x + y^2 + 3
    ω = [K(2), K(3)]
    @test [K(5), K(14), K(88)] == ESI.fast_multivariate_evaluate(R, f, ω, 3)
    @test [K(5), K(14)] == ESI.fast_multivariate_evaluate(R, f, ω, 2)
    @test [K(5)] == ESI.fast_multivariate_evaluate(R, f, ω, 1)
    @test [K(5), K(14), K(88), K(740), K(6580)] == ESI.fast_multivariate_evaluate(R, f, ω, 5)
    for K in [
            GF(2^31-1), 
            GF(4611686018427388039)
            # FqNmodFiniteField(fmpz(7), 5, :z)
        ]
        R, (x,y,z) = PolynomialRing(K, ["x","y","z"])
        cases = [
            (poly=x*y*z + y + z, ω=[K(2),K(3),K(5)]),
            (poly=x + y + z, ω=[K(2),K(3),K(5)])
        ]
        for case in cases
            poly, ω = case.poly, case.ω
            t = length(poly)
            T = rand((div(t, 2), t, rand(1:5)*t))
            evals_1 = ESI.fast_multivariate_evaluate(R, poly, ω, T)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:T-1)
            @test evals_1 == evals_true
        end
        ω = [K(5), K(3), K(2)]
        for i in 1:10
            poly = rand(R, 1:1000, 1:10)
            t = length(poly)
            T = rand((div(t, rand(2:5)), t, rand(1:5)*t))
            T == 0 && continue
            evals_1 = ESI.fast_multivariate_evaluate(R, poly, ω, T)
            evals_true = map(i -> evaluate(poly, ω .^ i), 0:T-1)
            @test evals_1 == evals_true
        end
    end

end
