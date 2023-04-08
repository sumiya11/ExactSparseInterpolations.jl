
include("../../src/ExactSparseInterpolations.jl")

using Nemo, Compose

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])
Ru, (x, u) = PolynomialRing(K, ["x","u"])

f = (x1 + x2^7 + x3^2*x4 + x4 + 8)*(x1 + x2^4 + x3*x4 + 3)*(x2^10 + x3^2 + 9)
@info "Degree in x2: $(degree(f, x2)), true number of factors: $(length(factor(f)))"
for _ in 1:100
    p = rand(K)
    f0 = evaluate(f, [R(p),x2,R(p),R(p)])
    fi = factor(f0)
    @info "Factored at $p: $(length(fi)) factors"
end

factor(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    skipcontent=true,
    benchmark=true,
    mainvar=:largedegree
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_main_var]
data = ExactSparseInterpolations.dump_data()
ExactSparseInterpolations.draw_graph(
    PDF("tree2.pdf", 16cm, 16cm), benchs[:v_tree], data...)
