ESI = ExactSparseInterpolations

@testset "generic" begin
    R, (x,) = PolynomialRing(QQ, ["x"])
    info = ESI.getboundsinfo(x^3 + 1)
    @test info.totaldeg == 3
    @test info.nterms == 2
    @test info.partialdegs == [3]

    R, (x,y,z) = PolynomialRing(GF(17), ["x","y","z"])
    info = ESI.getboundsinfo((x*y^2 - z)//(y^5 + y^2 + y + 1))
    @test info.numtotaldeg == 3
    @test info.numnterms == 2
    @test info.numpartialdegs == [1, 2, 1]
    @test info.dentotaldeg == 5
    @test info.dennterms == 4
    @test info.denpartialdegs == [0, 5, 0]
end
