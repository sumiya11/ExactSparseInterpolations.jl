using Nemo

# > make degree in some variable even smaller?

#=
    #1 Benchmark, in general

    Nice polynomial, 
    a couple of irreducible factors, 
    main variable of low degree,
    no transformation is needed
=#

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4,x5) = PolynomialRing(K, ["x$i" for i in 1:5])

f = (x1 + (x2 + x3 + x4 + x5 + 1)^5)*(x1^3 + (x2*x3 + x3*x4 + x2*x4 + x1*x5)^3 + 2);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_points_used]
benchs[:v_transform_degrees]
benchs[:t_total]
sum(benchs[:t_select_main_variable])
sum(benchs[:t_interpolation])
sum(benchs[:t_many_hensel_liftings])
sum(benchs[:t_evaluating_coefficients])
sum(benchs[:t_first_bivariate_factor])

##########################

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4,x5,x6,x7,x8) = PolynomialRing(K, ["x$i" for i in 1:8])

f = (x1 + (x2 + x3 + x4 + x5 + x6 + x7 + x8 + 1)^4)*(x1^2 + (x2*x3 + x3*x4 + x2*x4 + x2*x5 + x4*x7 + x5*x8 + 2)^2 + 2);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_points_used]
benchs[:v_transform_degrees]
benchs[:t_total]
sum(benchs[:t_select_main_variable])
sum(benchs[:t_interpolation])
sum(benchs[:t_many_hensel_liftings])
sum(benchs[:t_evaluating_coefficients])
sum(benchs[:t_first_bivariate_factor])

#=
    #2 Benchmark, searching for a power product

    Nice polynomial, 
    irreducible factors that depend on the same variable, 
    main variable of low degree
=#

# Trailing term of degree 6

# ssh demin@login.lix.polytechnique.fr
# ssh ron.medicis.polytechnique.fr
# TODO: Bug
K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4,x5) = PolynomialRing(K, ["x$i" for i in 1:5])

f = (x1 + (x2 + x3 + x4 + 1)^3)*(x1^5 + (x2*x3 + x3*x4 + x3*x5 + x2*x5 + x2*x4 + x4*x5)^2 + 1);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_points_used]
benchs[:v_transform_degrees]
benchs[:v_main_var]
benchs[:t_total]
sum(benchs[:t_select_main_variable])
sum(benchs[:t_interpolation])
sum(benchs[:t_many_hensel_liftings])
sum(benchs[:t_evaluating_coefficients])
sum(benchs[:t_first_bivariate_factor])

##########################
# Trailing term of degree 1

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4,x5) = PolynomialRing(K, ["x$i" for i in 1:5])

f = prod(x1^i + (x1 + x2 + x3 + x4 + x5 + (i - 1)) for i in 1:10);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_points_used]
benchs[:v_transform_degrees]
benchs[:v_main_var]
benchs[:t_total]
sum(benchs[:t_select_main_variable])
sum(benchs[:t_interpolation])
sum(benchs[:t_many_hensel_liftings])
sum(benchs[:t_evaluating_coefficients])
sum(benchs[:t_first_bivariate_factor])

##########################

#=
    #3 Benchmark, varying main var degree

    Nice polynomial, 
    irreducible factors that depend on the same variable, 
=#

K = Nemo.GF(2^62 + 135)
# K = Nemo.GF(4611686018427388319)
# K = Nemo.GF(fmpz(2)^120 + 451)
R, (x1,x2,x3) = PolynomialRing(K, ["x$i" for i in 1:3])

f = (x1^1 + (x2 + x3 + 1)^5)*(x1^(3) + (x2 + x3 + x2*x3 + 1)^6 + 2);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)

res = Dict();
benchs = ExactSparseInterpolations.dump_benchmarks()
for i in 1:14
    f = (x1^i + (x2 + x3 + 1)^5)*(x1^(i + 3) + (x2 + x3 + x2*x3 + 1)^6 + 2);
    length(f)

    fi = ExactSparseInterpolations.top_level_factorize(
        f, 
        benchmark=true,
        skipcontent=true
    )
    @assert prod(fi) == f "$i is bad"
    
    benchs = ExactSparseInterpolations.dump_benchmarks()

    if length(fi) > 2 || length(fi) == 1
        @warn "i = $i skipped"
        continue
    end
    @assert all(x->x.main_var_idx==1, benchs[:v_main_var])

    println("==============================")
    show(stdout, MIME"text/plain"(), benchs)
    println()

    res[3 + 2i] = (
        benchs[:t_total], 
        benchs[:t_evaluating_coefficients],
        benchs[:t_many_hensel_liftings], 
        benchs[:t_interpolation],
        benchs[:t_select_main_variable],
        benchs[:t_first_bivariate_factor] 
    )
end
res

using Plots
xs = sort(collect(keys(res)))
ys_total = [res[x][1][1] for x in xs]
ys_eval = [sum(res[x][2]) for x in xs]
ys_hensel = [sum(res[x][3]) for x in xs]
ys_interpolation = [sum(res[x][4]) for x in xs]
ys_mainvar = [sum(res[x][5]) for x in xs]
ys_first_bivariate_factor = [sum(res[x][6]) for x in xs]

begin
    pt = :tab10
    plo = plot(xs, ys_total, linewidth=2, label="total", palette=pt,
        xlabel="degree of main variable", ylabel="time (s)")
    plot!(xs, ys_mainvar, linewidth=2, label="find main variable", palette=pt)
    plot!(xs, ys_first_bivariate_factor, linewidth=2, label="first bivariate factoring", palette=pt)
    plot!(xs, ys_eval, linewidth=2, label="evaluating", palette=pt)
    plot!(xs, ys_hensel, linewidth=2, label="hensel lift", palette=pt)
    plot!(xs, ys_interpolation, linewidth=2, label="interpolation", palette=pt)
    savefig(plo, "benchmark_3.pdf")
end

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true
)

benchs = ExactSparseInterpolations.dump_benchmarks()

# ExactSparseInterpolations.select_main_variable(
#     (x1 - 1)*(x1 + x7 + 3)*(x1 - x6 + x7 + 2)
# )

# f = (x1 - 1)*(x1 + x7 + 3)*(x1 - x6 + x7 + 2)
f = (x1 + x2*x3 + x3*x4)*(x1 + 1)
fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true
)
ExactSparseInterpolations.top_level_factorize(f, benchmark=true)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_points_used]
benchs[:v_transform_degrees]
benchs[:v_transform_matrices]
benchs[:v_main_var]

factor(f)

##########################

K = Nemo.GF(2^31-1)
R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])
f = (x1 + (x2*x3 + x3*x4 + x4^2)^1)*((x1*(x3 + 1) + x2 + x3 + x4)^1 + 1);
fi = ExactSparseInterpolations.top_level_factorize(f, benchmark=true)
prod(fi) == f

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:t_total]
benchs[:v_main_var]

@time factor(f);

factor(f);
success, fi = ExactSparseInterpolations.find_some_factors(f, 100);
prod(fi) == f

ExactSparseInterpolations.factorize(f)

f = (x1 + (x2*x3 + x3*x4 + x4^2)^3)*(x1 + x4 + 1)*(x1 + (x2 + x3 + x4 + 1)^3)
factor(f)

ExactSparseInterpolations.factorize(f)

# f = prod(x1 + x2^rand(0:2) + x3^rand(0:2) + rand(K) for _ in 1:3)*prod(x2^2 + x3^rand(0:2) + x4^rand(0:2) + rand(K) for _ in 1:3)*prod(x4 + x5^rand(0:2) + x6^rand(0:2) + x7^rand(0:2) + rand(K) for _ in 1:3)

bench = false
@macroexpand ExactSparseInterpolations.@savetime bench :t_owo zz = 1 + 3
@macroexpand ExactSparseInterpolations.@savetime Main.bench :t_owo zz = 1 + 3
ExactSparseInterpolations.@savetime Main.bench :t_owo zz = 1 + 3
zz
ExactSparseInterpolations._factorize_benchmarks

##### 

K = Nemo.GF(2^31-1)
R, (x1,x2,x3) = PolynomialRing(K, ["x$i" for i in 1:3])
f = (x1 + (x2*x3 + x3^2)^2 + 99)*((x1*(x3 + 1) + x2 + x3)^1 + 1);
g = (x1 + (x2^3*x3 + x3)^3 + 5)*((x1 + 5x3 + 8)^1);
h = (x1 + (x2 + x3 + 2)^3 + 88)*(x1 + x2 + x3 + 3);
f = f*g*h;
fi = ExactSparseInterpolations.top_level_factorize(f, benchmark=true);
prod(fi) == f

# not squarefree
factor(evaluate(f, [x1, R(0), R(0)]))

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:t_total]
benchs[:t_many_hensel_liftings]
sum(benchs[:t_removing_content])

length(factor(f)), length(fi)
