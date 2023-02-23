ESI = ExactSparseInterpolations

@testset "hensel lifting Zp[x] -> Zp[x,y]" begin
    R, (x, y) = PolynomialRing(GF(2^31-1), ["x", "y"])
    
    cases = [
        (m=y, l=2, Fi = [x + y, x^2 + 1], lc=R(1)),
        (m=y, l=2, Fi = [x + y, x^2 + 2*y + 8], lc=R(1)),
        (m=y, l=3, Fi = [x + y, x^2 + 2*y + 1, x + 3*y^2 + 2], lc=R(1)),
        (m=y, l=2, Fi = [x + y, x + 2], lc=y + 1),
        (m=y, l=15, Fi = [x + 1, x + 100y^10], lc=y^4 + 1),
    ]

    for case in cases
        m, l, Fi, lc = case.m, case.l, case.Fi, case.lc
        fi = map(f -> mod(f, m), Fi)
        F = lc*prod(Fi)
        fi_star = ESI.hensel_multifactor_lifting(F, fi, l, m)
        @test all(f_star -> isone(ESI.leading_coefficient_in(f_star, x)), fi_star)
        @test all(f1f2 -> mod(f1f2[1], m) == mod(f1f2[2], m), zip(fi, fi_star))
        A = ESI.leading_coefficient_in(F, x)*prod(fi_star)
        @test mod(lc*A, m^l) == mod(lc*F, m^l)
        @test all(in(fi_star), Fi) && all(in(Fi), fi_star)
    end
end
