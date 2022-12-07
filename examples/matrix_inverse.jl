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

function direct_inverse(R, M)
    inv(M)
end

function interpolated_inverse(R, M, Minv)
    evaluationmap = point -> inv(map(ij -> evaluate(ij, point), M))
    cij = [
        [point -> evaluationmap(point)[i, j] 
        for i in 1:size(M, 1)]
            for j in 1:size(M, 2)
    ]
    Iij = [
        [ExactSparseInterpolations.vanDerHoevenLecerf(R, maxdeg(Minv[i, j])...)
        for i in 1:size(M, 1)]
            for j in 1:size(M, 2)
    ]
    ans = zero(M)
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            P, Q = ExactSparseInterpolations.interpolate!(Iij[i][j], cij[i][j])
            ans[j, i] = P // Q
        end
    end
    ans
end

function benchmark(R, ns; 
        term_range=0:2, 
        exponent_range=0:2, 
        coeff_range=1:10)
    params = (term_range, exponent_range, coeff_range)
    for n in ns
        while true
            M = generate_matrix(R, n, params...)
            iszero(det(M)) && continue
            Minv1 = direct_inverse(R, M)
            Minv2 = interpolated_inverse(R, M, Minv1)
            # @info "" Minv1
            # @info "" Minv2
            @assert Minv1 == Minv2
            b1 = @benchmarkable direct_inverse($R, $M) samples=3
            b2 = @benchmarkable interpolated_inverse($R, $M, $Minv1) samples=3
            println(run(b1))
            println(run(b2))
            println("-----------")
            break
        end
    end
end

R, (x, y) = PolynomialRing(QQ, ["x", "y"])

a = benchmark(R, 2:8)
