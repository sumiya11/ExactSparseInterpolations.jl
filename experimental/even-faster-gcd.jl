using Nemo
using BenchmarkTools

# from nanoseconds to seconds
sround(t) = t # round(t, digits=7);
ns2s(t) = sround(t / 1e9);

function benchmark_gcd(GCD, ns)
    R, x = PolynomialRing(GF(2^31-1), "x")
    times = Float64[]
    for i in ns
        n = 2^i
        a, b, c = div(n, 3), div(n, 6), div(n, 2)
        @assert a + b + c <= n
        # f = (x - 1)^a*(x + 4)^b*(x + 3)^c
        # g = (x - 1)^b*(x + 3)^a*(x + 7)^c
        f = (x + 1)^n
        g = (x + 2)^(n - 1)
        bench = @benchmarkable $GCD($f, $g) samples=3
        push!(times, ns2s(minimum(run(bench)).time))
    end
    times
end

ns = 7:18
esi_fast_gcd = benchmark_gcd(ExactSparseInterpolations.fastgcd, ns);
@info "" esi_fast_gcd

# esi_slow_gcd = benchmark_gcd(ExactSparseInterpolations.slowgcd, ns);
# @info "" esi_slow_gcd

## For threshold = 0
# 12-element Vector{Float64}:
# │      0.0009267
# │      0.0018341
# │      0.0040236
# │      0.0094329
# │      0.0222495
# │      0.0516653
# │      0.1278791
# │      0.3335785
# │      0.85574
# │      2.0189089
# │      4.6595947
# └     10.6043042

## For threshold = 2^2
# │    12-element Vector{Float64}:
# │      0.0006336
# │      0.0015129
# │      0.0032953
# │      0.0077897
# │      0.0189591
# │      0.0458171
# │      0.1156284
# │      0.2751452
# │      0.7140653
# │      1.9786393
# │      4.5539313
# └     10.2641142

## For threshold = 2^4
# 12-element Vector{Float64}:
# │     0.0004498
# │     0.0010499
# │     0.0024996
# │     0.0070707
# │     0.0182644
# │     0.0389071
# │     0.0979718
# │     0.2468699
# │     0.6161165
# │     1.4331917
# │     3.6708854
# └     8.4749868

