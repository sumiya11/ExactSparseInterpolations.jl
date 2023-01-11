using Nemo, BenchmarkTools
# using ExactSparseInterpolations

### Over Z/Zp ###

F = GF(2^16 + 1)
R, (x,) = PolynomialRing(F, ["x",])

poly = (x^10 + 5x^9 + 11x + 18)^100;
@show length(poly)
f = ExactSparseInterpolations.Blackbox(poly);

bot = ExactSparseInterpolations.FasterBenOrTiwari(R, length(poly))
@time interpolated = ExactSparseInterpolations.interpolate!(bot, f);

@assert interpolated == poly
