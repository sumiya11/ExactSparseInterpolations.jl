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

f = (x1 + (x2 + x3 + x4 + x5 + 1)^2)*(x1^3 + (x2*x3 + x3*x4 + x2*x4 + x1*x5)^1 + 2);
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
# Trailing term of degree 6

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4,x5) = PolynomialRing(K, ["x$i" for i in 1:5])

f = (x1 + (x2 + x3 + x4 + 1))*(x1 + (x2*x3 + x3*x4 + x3*x5 + x2*x5 + x2*x4 + x4*x5)^2);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_main_var]
benchs[:v_transform_degrees]
benchs[:v_transform_trailing_term]
benchs[:v_transform_matrices]

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

f = (x1^13 + (x2 + x3 + 1)^5)*(x1^(3 + 13) + (x2 + x3 + x2*x3 + 1)^6 + 2);
length(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)

Ru,(t, u) = K["t","u"]
F_sub_u = evaluate(f, [t, u, u])
@time ExactSparseInterpolations.revealing_bivariate_factorization_ff(F_sub_u);

factor(evaluate(F, [t,Ru(0)]))

F = F_sub_u
@time begin
    R = parent(F)
    K, (t, u) = base_ring(R), gens(R)
    n, d = degree(F, t), degree(F, u)
    b = ExactSparseInterpolations.leading_coefficient_in(F, t)
    # modulo for reduction to univariate case
    m = u
    f = mod(F, m)
    @assert n >= 1
    @assert isone(gcd(f, derivative(f, t)))
    @assert !iszero(mod(b, m))
    Runiv, _ = K["t"] 
    # l is the bound on the length of hensel iteration
    l = d + 1 + degree(b, u)
    # factor F mod m
    # into factors f1..fr
    funiv = to_univariate(Runiv, f)
    factorization = Nemo.factor(funiv)
    fiuniv = collect(keys(factorization.fac))
    # lift the factors f1..fr via Hensel lifting
    fi = map(f -> evaluate(f, t), fiuniv)
    filifted = ExactSparseInterpolations.hensel_multifactor_lifting(F, fi, l, m)

    Fi, Si = ExactSparseInterpolations.recombine_factors(F, filifted, m^l)
end;

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

xs0 = deepcopy(xs)
res0 = deepcopy(res)

xs = xs[1:end-2]
ys_total = ys_total[1:end-2]
ys_eval = ys_eval[1:end-2]
ys_hensel = ys_hensel[1:end-2]
ys_interpolation = ys_interpolation[1:end-2]
ys_mainvar = ys_mainvar[1:end-2]
ys_first_bivariate_factor = ys_first_bivariate_factor[1:end-2]
begin
    pt = :tab10
    plo = plot(xs, ys_total, linewidth=2, label="total", palette=pt,
        xlabel="degree of main variable", ylabel="time (s)")
    plot!(xs, ys_mainvar, linewidth=2, label="find main variable", palette=pt)
    plot!(xs, ys_first_bivariate_factor, linewidth=2, label="first bivariate factoring", palette=pt)
    plot!(xs, ys_eval, linewidth=2, label="evaluating", palette=pt)
    plot!(xs, ys_hensel, linewidth=2, label="hensel lift", palette=pt)
    plot!(xs, ys_interpolation, linewidth=2, label="interpolation", palette=pt)
    savefig(plo, "benchmark_4.pdf")
end


##########################

#=
    #4 Benchmark, varying the number of terms
=#

using Nemo
K = Nemo.GF(2^62 + 135)
# K = Nemo.GF(4611686018427388319)
# K = Nemo.GF(fmpz(2)^120 + 451)
R, xs = PolynomialRing(K, ["x$i" for i in 1:8], ordering=:degrevlex)

xso = xs[2:end]
f = (xs[1] + (sum(xso) + 1)^3)*(xs[1] + (sum(xso)^4));
length(f)

@profview fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
);
length(fi)
prod(fi) == f

res = Dict();
ts = []
benchs = ExactSparseInterpolations.dump_benchmarks()
for i in 2:11
    K = Nemo.GF(2^62 + 135)
    R, xs = PolynomialRing(K, ["x$i" for i in 1:8])

    d1 = div(i, 2)
    d2 = i - d1
    xso = xs[2:end]
    f = (xs[1] + (sum(xso) + 1)^d1)*(xs[1] + (sum(xso)^d2));
    t = length(f)
    push!(ts, t)

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

    println("==============================")
    show(stdout, MIME"text/plain"(), benchs)
    println()

    res[i] = (
        benchs[:t_total],
        benchs[:t_evaluating_coefficients],
        benchs[:t_many_hensel_liftings], 
        benchs[:t_interpolation],
        benchs[:t_select_main_variable],
        benchs[:t_first_bivariate_factor],
        benchs[:t_postprocessing],
        benchs[:t_exact_division],
        benchs[:t_find_power_product],
    )
end
res

using Plots
begin
    xs = sort(collect(keys(res)))
    ys_total = [sum(res[x][1]) for x in xs]
    ys_eval = [sum(res[x][2]) for x in xs]
    ys_hensel = [sum(res[x][3]) for x in xs]
    ys_interpolation = [sum(res[x][4]) for x in xs]
    ys_mainvar = [sum(res[x][5]) for x in xs]
    ys_first_bivariate_factor = [sum(res[x][6]) for x in xs]
    ys_postprocessing = [sum(res[x][7]) for x in xs]
    ys_exact_division = [sum(res[x][8]) for x in xs]
    ys_find_power_product = [sum(res[x][9]) for x in xs]
    ys_almost_total = ys_eval .+ ys_hensel .+ ys_interpolation .+ ys_mainvar .+ ys_first_bivariate_factor .+ ys_postprocessing .+ ys_exact_division .+ ys_find_power_product
end

begin
    prcnt = round(Int, 100*sum(ys_almost_total) / sum(ys_total))
    title = "8 variables, recorded $prcnt% of runtime"
    pt = :tab10
    # plo = plot(xs, ys_total, linewidth=2, label="total", palette=pt,
    #     xlabel="# terms in input", ylabel="time (s)", title=title,
    #     xticks=(1:length(ts), ts))
    plot(xs, ys_find_power_product, linewidth=2, 
        label="find power product", palette=pt,
        xlabel="# terms in input", ylabel="time (s)", title=title,
        xticks=(1:length(ts), ts))
    plot!(xs, ys_mainvar, linewidth=2, label="find main variable", palette=pt)
    plot!(xs, ys_first_bivariate_factor, linewidth=2, label="first bivariate factoring", palette=pt)
    plot!(xs, ys_eval, linewidth=2, label="evaluation", palette=pt)
    plot!(xs, ys_hensel, linewidth=2, label="hensel lift", palette=pt)
    plot!(xs, ys_interpolation, linewidth=2, label="interpolation", palette=pt)
    plot!(xs, ys_postprocessing, linewidth=2, label="postprocessing", palette=pt)
    plot!(xs, ys_exact_division, linewidth=2, label="exact division", palette=pt)
    # savefig(plo, "benchmark_4.pdf")
end


#=
    #4 Benchmark, exhaustive search for a better matrix
=#

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])

f = (x1 + (x2^2 + x3^2 + x4^2)^3)*(x1 + (x2 + x3 + x4 + 1)^3);
monoms_no_x1 = filter(m -> degree(m, x1) == 0, collect(monomials(f)))
length(f), minimum(total_degree, monoms_no_x1)
filter(m -> total_degree(m) == minimum(total_degree, monoms_no_x1), (collect(monomials(f))))
f;

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true,
    monomtransform=:exhaustive
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_main_var]
benchs[:v_transform_degrees]
benchs[:v_transform_trailing_term]
benchs[:v_transform_matrices]
benchs[:v_transform_matrices][1][2]
benchs[:v_transform_matrices][2][2]
benchs[:v_transform_matrices][3][2]


K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])

f = (x1 + (x2^2 + x3 + x4)^5)*(x1 + (x2 + x3 + x4 + 1)^5);
monoms_no_x1 = filter(m -> degree(m, x1) == 0, collect(monomials(f)))
length(f), minimum(total_degree, monoms_no_x1)
filter(m -> total_degree(m) == minimum(total_degree, monoms_no_x1), (collect(monomials(f))))
f;

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_main_var]
benchs[:v_transform_degrees]
benchs[:v_transform_trailing_term]
benchs[:v_transform_matrices]

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

#=
    #5 Benchmark.
    Using many points for bivariate factorization
    + check if small degrees are better than large ones.
=#

using Nemo

K = Nemo.GF(2^62 + 135)
R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])
Ru, (x, u) = PolynomialRing(K, ["x","u"])

f = (x1^2 + x2*x3 + x3*x4 - 1)*(x1 + 2x2)*(x1 + x2^3*x3^4 + 8)*(x1 + (x3 + 9)^2)
factor(f)

fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    skipcontent=true,
    benchmark=true,
)
prod(fi) == f
map(length, fi)

data = ExactSparseInterpolations.dump_data()
benchs = ExactSparseInterpolations.dump_benchmarks()
ExactSparseInterpolations.draw_graph(benchs[:v_tree], data...)

benchs[:v_main_var]
benchs[:v_tree]

f = (x1 + x2 + x3)*(x1^3 + x3)*(x1*x2 + x3 + x4)
fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    skipcontent=true,
    benchmark=true,
    beautifuly=true
)
prod(fi) == f
data = ExactSparseInterpolations.dump_data()
benchs = ExactSparseInterpolations.dump_benchmarks()
ExactSparseInterpolations.draw_graph(benchs[:v_tree], data...)


Fi

F_u_c2 = evaluate(f, [x, 2u, 3u, 5u])
factor(F_u_c2)

p = K(0)
fi = map(f -> evaluate(f, [x, p]), Fi)
m = u - p
l = degree(F_u_c2, u) + 1
filifted = ExactSparseInterpolations.hensel_multifactor_lifting(F_u_c2, fi, l, m)



f = (x1^15 + (x2*x3^5 + x2*x3^5 + x2*x4^5 + 2))*(x1^13 + (x2^3 + x3^4 + x4^2 + 1));

fu = evaluate(f, [x1,x2,x2,x2])
factor(fu)

fu_0 = evaluate(f, [x1,K(0),K(0),K(0)])
length(factor(fu_0))
fu_1 = evaluate(f, [x1,K(1),K(1),K(1)])
length(factor(fu_1))
fu_235 = evaluate(f, [x1,K(2),K(3),K(5)])
length(factor(fu_235))
fu_5 = evaluate(f, [x1,K(5),K(5),K(5)])
factor(fu_5)

for i in 1:100
    a,b,c = rand(K),rand(K),rand(K)
    a = b = c
    fu_abc = evaluate(f, [x1,a*rand(K),b*rand(K),c*rand(K)])
    @info length(factor(fu_abc))
end

# Substitute c ^ 1
fu_0 = evaluate(f, [x1,x2,x2,x2])
factor(fu_0)

@profview fi = ExactSparseInterpolations.top_level_factorize(
    f, 
    benchmark=true,
    skipcontent=true
)
prod(fi) == f
map(length, fi)

benchs = ExactSparseInterpolations.dump_benchmarks()
benchs[:v_main_var]
benchs[:v_transform_degrees]
benchs[:v_transform_trailing_term]
benchs[:v_transform_matrices]
benchs[:t_first_bivariate_factor]

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
