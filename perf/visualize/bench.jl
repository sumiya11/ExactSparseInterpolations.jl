using BenchmarkTools, Random
using Nemo, Primes
using ExactSparseInterpolations
# for saving results into a file
using JLD

include("../testsuite.jl")

BenchmarkTools.DEFAULT_PARAMETERS.samples = 3
Random.seed!(888)

# Benchmark the interpolation of random functions over a finite field
# of small characteristic
big_prime = 2^62 + 135
@assert Primes.isprime(big_prime)
ground_field = GF(big_prime)

function benchmark_function_vdhl(func)
    R = base_ring(parent(func))
    metainfo = ExactSparseInterpolations.getboundsinfo(func)
    blackbox = ExactSparseInterpolations.Blackbox(func)
    vdhl = ExactSparseInterpolations.FasterVanDerHoevenLecerf(R, metainfo)    
    P, Q = ExactSparseInterpolations.interpolate!(vdhl, blackbox)
    @assert P//Q == func
    runtime = @belapsed ExactSparseInterpolations.interpolate!($vdhl, $blackbox) setup=(blackbox = ExactSparseInterpolations.Blackbox($func))
    runtime
end

@info "Benchmarking random functions over $ground_field"
nterms = [floor(Int, 1.5^i) for i in 5:18]  # from 7 to 1.5k
nvars = 2:1:8
degrees = [8, 16, 32, 64]
coeff_range = 1:2^16
data = Array{Float64, 3}(undef, length(nterms), length(nvars), length(degrees))

@info "Terms: $nterms, vars: $nvars, degrees: $degree, coeff_range: $coeff_range"
for (i, nt) in enumerate(nterms)
    for (j, nv) in enumerate(nvars)
        for (k, nd) in enumerate(degrees)
            num_degree = den_degree = div(nd, 2)
            num_terms = den_terms = div(nt, 2)
            func = random_rational_function(
                coeff_range, nv, num_degree, den_degree, 
                num_terms, den_terms, 
                ground_field=ground_field
            )
            @info "Function ($nt, $nv, $nd): $(ExactSparseInterpolations.getboundsinfo(func))"
            runtime = benchmark_function_vdhl(func)
            @info "Runtime: $runtime"
            data[i, j, k] = runtime
        end
    end
end

# save the tensor with results
save((@__DIR__) * "/data.jld", "data", data)
