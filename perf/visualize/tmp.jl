using BenchmarkTools, Random
using Nemo, Primes
using ExactSparseInterpolations

include("../testsuite.jl")

F = GF(2^30*3 + 1)

nterms = 100
nv = 5
deg = 16

func = random_rational_function(1:2^16, nv, deg, deg, nterms, nterms, ground_field=F);
R = base_ring(parent(func))
metainfo = ExactSparseInterpolations.getboundsinfo(func)
blackbox = ExactSparseInterpolations.Blackbox(func);

vdhl = ExactSparseInterpolations.FasterVanDerHoevenLecerf(R, metainfo);

# @profview ExactSparseInterpolations.interpolate!(vdhl, blackbox);

# Given the function `GCD` and the range `ns`
# benchmark the function `GCD` on polynomials of degrees 2^n for n in `ns` 
# and return the array of running times
function benchmark_gcd(GCD, ns)
    R, x = PolynomialRing(GF(2^31-1), "x")
    times = Float64[]
    for i in ns
        n = 2^i
        a, b, c = div(n, 3), div(n, 6), div(n, 2)
        @assert a + b + c <= n
        f = (x - 1)^a*(x + 4)^b*(x + 3)^c
        g = (x - 1)^b*(x + 3)^a*(x + 7)^c
        bench = @benchmarkable $GCD($f, $g) samples=3
        push!(times, minimum(run(bench)).time)
    end
    times
end
esi_fast_gcd = benchmark_gcd(ExactSparseInterpolations.fastgcd, ns);
