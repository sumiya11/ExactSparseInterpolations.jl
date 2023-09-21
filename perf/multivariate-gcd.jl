using ExactSparseInterpolations, Nemo
using BenchmarkTools

R, x = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:10]...])

P = (x[1] + x[2] + 2) * sum(x)^2 * (prod(x) - 1)
Q = (x[1] + x[2] + 3) * sum(x)^3 * (prod(x) + 1)

G1 = Nemo.gcd(P, Q)
G2 = multivariate_gcd(P, Q)
