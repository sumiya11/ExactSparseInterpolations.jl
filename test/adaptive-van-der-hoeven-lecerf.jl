ESI = ExactSparseInterpolations

maxdeg(f) = maximum(map(t -> total_degree(t), terms(f)), init=0) 

@testset "simultaneous adaptive v-d-h-l" begin
    R, (x1,x2,x3) = PolynomialRing(GF(fmpz(4611686018427388039)), ["x1","x2","x3"])
    case = [
        [R(1)//R(x2^3), R(1)//(x2), (x1*x2*x3)//(x2 + 8)],
        [R(1)//R(1), R(1)//(x2^4), (x1 - x2)^2 // (x2 - x3)^2],
        [R(4)//R(5), x2^4//x1^3]
    ]
    fs = [
        [ESI.Blackbox(c) for c in casei]
        for casei in case
    ]
    ans = [
        [one(R)//one(R) for c in casei]
        for casei in case
    ]
    Is = ESI.several_adaptiveVanDerHoevenLecerf(R, 4, 4, fs)
    i = 0
    all_interpolated = false
    while !all_interpolated
        all_interpolated = true
        x = ESI.next_point!(Is)
        for i in 1:length(Is)
            for j in 1:length(Is[i])
                flag, (P, Q) = ESI.next!(Is[i][j], fs[i][j](x))
                all_interpolated = all_interpolated && flag
                ans[i][j] = P//Q
            end
        end
    end
    @test case == ans

    case = [
        [(x1 - x2)^8//x2^3, R(1)//R(99)],
        [R(999)//(x1*x2*x3)^3, (x1 + x2 + x3 + 1)^3//R(8)],
    ]
    fs = [
        [ESI.Blackbox(c) for c in casei]
        for casei in case
    ]
    ans = [
        [one(R)//one(R) for c in casei]
        for casei in case
    ]
    Is = ESI.several_adaptiveVanDerHoevenLecerf(R, 8, 9, fs)
    i = 0
    all_interpolated = false
    while !all_interpolated
        all_interpolated = true
        x = ESI.next_point!(Is)
        for i in 1:length(Is)
            for j in 1:length(Is[i])
                flag, (P, Q) = ESI.next!(Is[i][j], fs[i][j](x))
                all_interpolated = all_interpolated && flag
                ans[i][j] = P//Q
            end
        end
    end
    @test case == ans
end

@testset "adaptive van-der-Hoeven-Lecerf, Q" begin
    R, (x1,x2) = PolynomialRing(QQ, ["x1","x2"])
    cases = [
        (x1)//(2x2 + 1), R(2)//R(3), R(0)//R(1),
        (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3,
        (x1 + x2)^5//(x1 - x2)^6, (x1*x2 - x1 - x2)^3 // R(8)
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        t1, t2 = length(numerator(case)), length(denominator(case))
        vdhl = ESI.adaptiveVanDerHoevenLecerf(R, n, d)
        flag = false
        P, Q = one(R), one(R)
        for i in 1:(2*max(t1, t2)+1)*(n + d + 2)
            x = ESI.next_point!(vdhl)
            y = f(x)
            flag, (P, Q) = ESI.next!(vdhl, y)
        end
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
        for i in 1:(n + d + 2)
            x = ESI.next_point!(vdhl)
            y = f(x)
            flag, (P, Q) = ESI.next!(vdhl, y)
        end
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

    f = ESI.Blackbox((x1 + x2 - x1*x2)//(x1 + x2))
    vdhl = ESI.vanDerHoevenLecerf(R, 3, 2)
    (P, Q) = ESI.interpolate!(vdhl, f)
    @test P//Q == (x1 + x2 - x1*x2)//(x1 + x2)

    R, xi = PolynomialRing(QQ, ["x$i" for i in 1:10])
    cases = [
        sum(xi) // (prod(xi)),
        (xi[1] + xi[10])//(xi[1] - xi[10] - 2),
        (xi[1] + xi[3] + xi[5] + xi[7])^5//(xi[1] + 2),
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        t1, t2 = length(numerator(case)), length(denominator(case))
        vdhl = ESI.adaptiveVanDerHoevenLecerf(R, n, d)
        flag = false
        P, Q = one(R), one(R)
        for i in 1:(2*max(t1, t2)+1)*(n + d + 2)
            x = ESI.next_point!(vdhl)
            y = f(x)
            flag, (P, Q) = ESI.next!(vdhl, y)
        end
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

end

@testset "adaptive van-der-Hoeven-Lecerf, Z/Zp" begin
    R, (x1,x2) = PolynomialRing(GF(2^31-1), ["x1","x2"])
    cases = [
        (x1)//(2x2 + 1), R(2)//R(3), R(0)//R(1),
        (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        t1, t2 = length(numerator(case)), length(denominator(case))
        vdhl = ESI.adaptiveVanDerHoevenLecerf(R, n, d)
        flag = false
        P, Q = one(R), one(R)
        for i in 1:(2*max(t1, t2)+1)*(n + d + 2)
            x = ESI.next_point!(vdhl)
            y = f(x)
            flag, (P, Q) = ESI.next!(vdhl, y)
        end
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
        for i in 1:(n + d + 2)
            x = ESI.next_point!(vdhl)
            y = f(x)
            flag, (P, Q) = ESI.next!(vdhl, y)
        end
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

    R, xi = PolynomialRing(GF(fmpz(2)^300 + 157), ["x$i" for i in 1:10])
    cases = [
        (xi[1] - xi[2])^40//(2xi[1]*xi[2]),
        sum(xi) // (prod(xi)),
        (xi[1] + xi[10])//(xi[1] - xi[10] - 2),
        (xi[1] + xi[3] + xi[5] + xi[7])^5//(xi[1] + 2),
    ]
    for case in cases
        f = ESI.Blackbox(case)
        n, d = maxdeg(numerator(case)), maxdeg(denominator(case))
        t1, t2 = length(numerator(case)), length(denominator(case))
        vdhl = ESI.adaptiveVanDerHoevenLecerf(R, n, d)
        flag = false
        P, Q = one(R), one(R)
        for i in 1:(2*max(t1, t2)+1)*(n + d + 2)
            x = ESI.next_point!(vdhl)
            y = f(x)
            flag, (P, Q) = ESI.next!(vdhl, y)
        end
        @test isone(trailing_coefficient(Q))
        @test P//Q == case
    end

end
