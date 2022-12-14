using Nemo
using ExactSparseInterpolations

### Over F_2^16 ###

F = FqNmodFiniteField(fmpz(2), 16, :z)
R, (x,) = PolynomialRing(F, ["x",])

poly = (x^10 + x^9 + x + 1)^100
f = ExactSparseInterpolations.Blackbox(poly)

bot = ExactSparseInterpolations.FasterBenOrTiwari(R, length(poly))
@time interpolated = ExactSparseInterpolations.interpolate!(bot, f)
#  0.087853 seconds (49.17 k allocations: 1.548 MiB)

@assert interpolated == poly

### Over Z/Zp ###

F = GF(2^16 + 1)
R, (x,) = PolynomialRing(F, ["x",])

poly = (x^10 + 5x^9 + 11x + 18)^100
f = ExactSparseInterpolations.Blackbox(poly)

bot = ExactSparseInterpolations.FasterBenOrTiwari(R, length(poly))
@time interpolated = ExactSparseInterpolations.interpolate!(bot, f);
#  0.452306 seconds (1.95 M allocations: 36.390 MiB, 13.25% gc time)

@assert interpolated == poly
