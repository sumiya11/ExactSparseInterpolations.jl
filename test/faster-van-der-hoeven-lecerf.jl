ESI = ExactSparseInterpolations

@testset "fast vdH-Lecerf, Zp and Fq^n" begin
    for ground in [
        Nemo.GF(3*2^30+1), 
        Nemo.FqNmodFiniteField(fmpz(5), 12, :z)
    ]
        R, (x1,x2) = PolynomialRing(ground, ["x1","x2"])
        cases = [
            1//(x1*x2), 1//(x1*x2)^3,
            (x1)//(2x2 + 1), R(2)//R(3),
            (3x1 + 4x2)//(x2), (x1 + x2)^3//(x1 - x2)^3,
            (x1 + x2)^4//(x1 - x2)^6, (x1*x2 - x1 - x2)^3 // R(3),
            (x1 + x2 + 1)^30//(x1 - x2),
            (x1^(2^10) + x2^(2^5) + 1)//R(8)
        ]
        for case in cases
            f = ESI.Blackbox(case)
            info = ESI.getboundsinfo(case)
            @info "" f info
            vdhl = ESI.FasterVanDerHoevenLecerf(R, info)
            P, Q = ESI.interpolate!(vdhl, f)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
        end

        R, xi = PolynomialRing(ground, ["x$i" for i in 1:10])
        cases = [
            sum(xi) // (prod(xi)),
            (xi[1] + xi[10])//(xi[1] - xi[10] - 2),
            (xi[1] + xi[3] + xi[5] + xi[7])^5//(xi[1] + 2),
        ]
        for case in cases
            f = ESI.Blackbox(case)
            info = ESI.getboundsinfo(case)
            vdhl = ESI.FasterVanDerHoevenLecerf(R, info)
            P, Q = ESI.interpolate!(vdhl, f)
            @test isone(trailing_coefficient(Q))
            @test P//Q == case
        end
    end
end

@testset "fast vdH-Lecerf, warns and errors" begin
    ESI.Random.seed!(1)

    ground = GF(101)
    R, (x,y,z) = PolynomialRing(ground, ["x","y","z"])
    cases = [
        (
            func=(x - y) // (y - z),
            ans=:works
        ),
        (
            # total degree in kronecker is 80, which is < 101
            func=(x*y*z^3) // (R(1)),
            ans=:works
        ),
        (
            # total degree in kronecker is 98, which is < 101
            func=(y^7 + z) // (R(1)),
            ans=:works
        ),
        (
            # total degree in kronecker is 120, which is > 120
            func=(x*y*z^4) // (R(1)),
            ans=:fails
        ),
        # (
        #     func=(x^10 + y^10 + z^10)//(x + y),
        #     ans=:warn
        # ),
        # (
        #     func=(x^2 + y^2 + z^2 + x*y + x*z + y*z + x + y + z + 1)//(x),
        #     ans=:err
        # )
    ]
    # for case in cases
    #     bb = ESI.Blackbox(case.func)
    #     info = ESI.getboundsinfo(case.func)
    #     vdhl = ESI.FasterVanDerHoevenLecerf(R, info)
    #     println(vdhl)
    #     if case.ans === :works
    #         P, Q = ESI.interpolate!(vdhl, bb)
    #         @test P//Q == case.func
    #     elseif case.ans === :fails
    #         P, Q = ESI.interpolate!(vdhl, bb)
    #         @test P//Q != case.func
    #     elseif case.ans === :warn
    #         @test_logs (:warn,) ESI.interpolate!(vdhl, bb)
    #     elseif ans === :err
    #         @test_throws AssertionError ESI.interpolate!(vdhl, bb)
    #     end
    # end
end
