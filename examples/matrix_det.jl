include("../src/ExactSparseInterpolations.jl")

using Nemo
using BenchmarkTools

maxdeg(f) = maximum(map(t -> total_degree(t), terms(f)), init=0) 
maxdeg(f::Nemo.Generic.Frac) = maxdeg(numerator(f)), maxdeg(denominator(f))

function generate_matrix(R, n, term_range, exponent_range, coeff_range)
    S = MatrixSpace(FractionField(R), n, n)
    A = zero(S)
    for i in 1:n
        for j in 1:n
            A[i, j] = rand(R, term_range, exponent_range, coeff_range)
        end
    end
    A
end

function direct_det(R, M)
    det(M)
end

function interpolated_det(R, M, d)
    c = point -> det(map(ij -> evaluate(ij, point), M))
    I = ExactSparseInterpolations.vanDerHoevenLecerf(R, maxdeg(d)...)
    P, Q = ExactSparseInterpolations.interpolate!(I, c)
    P // Q
end

function benchmark(R, ns; 
        term_range=0:2, 
        exponent_range=0:1, 
        coeff_range=1:1)
    params = (term_range, exponent_range, coeff_range)
    for n in ns
        while true
            M = generate_matrix(R, n, params...)
            iszero(det(M)) && continue
            d1 = direct_det(R, M)
            b1 = @benchmarkable direct_det($R, $M) samples=3
            b2 = @benchmarkable interpolated_det($R, $M, $d1) samples=3
            println(run(b1))
            println(run(b2))
            println("-----------")
            break
        end
    end
end

p = BigInt(2)^64 - 59
R, (x, y, z, t) = PolynomialRing(GF(fmpz(p)), ["x", "y","z","t"])

a = benchmark(R, 2:8)
