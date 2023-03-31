using Nemo

K = Nemo.GF(2^31-1)
R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])
xs = gens(R)

f = (x1 + 1)*(x2 + x3)

foo = evaluate(f, [x1, x1*x2, x3, x4])

fi = ExactSparseInterpolations.top_level_factorize(
    foo, 
    benchmark=true
)

factor(foo)

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
