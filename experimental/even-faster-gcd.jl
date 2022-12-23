using Nemo
using BenchmarkTools

# from nanoseconds to seconds
sround(t) = round(t, digits=4);
ns2s(t) = sround(t / 1e9);

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
        push!(times, ns2s(minimum(run(bench)).time))
    end
    times
end

ns = 7:18
esi_fast_gcd = benchmark_gcd(ExactSparseInterpolations.fastgcd, ns);
@info "" esi_fast_gcd
