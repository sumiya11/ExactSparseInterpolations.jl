using Nemo, Plots
using ProgressMeter
using Logging, Statistics, Combinatorics
using StatsBase, Measures
using LaTeXStrings

global_logger(ConsoleLogger(Logging.Warn))
const _progressbar_color = :magenta
const _progressbar_value_color = :magenta
progressbar_enabled() = Logging.Info <= Logging.min_enabled_level(current_logger())

function number_of_monoms_of_degree(n, d)
    binomial(BigInt(n + d - 1), BigInt(d))
end

function number_of_monoms_up_to_degree(n, d)
    binomial(BigInt(n + d), BigInt(d))
end

function random_poly_of_total_degree(R, T, D)
    nmons = [number_of_monoms_of_degree(nvars(R), i) for i in 0:D]
    weights = nmons ./ sum(nmons)
    poly = zero(R)
    K = base_ring(R)
    function rand_vector_with_sum(n, s)
        res = zeros(Int, n)
        while sum(res) < s
            res[rand(1:n)] += 1
        end
        res
    end
    for i in 1:T
        cf = rand(K)
        d = sample(0:D, Weights(weights))
        exp_vector = rand_vector_with_sum(Nemo.nvars(R), d)
        poly += R([cf], [exp_vector])
    end
    poly
end

function random_poly_of_density(R, D, density)
    T = 1
    f = R(rand(base_ring(R)))
    while length(f) / number_of_monoms_up_to_degree(nvars(R), D) < density
        T = 2T
        f = random_poly_of_total_degree(R, T, D)
    end
    f
end

function benchmark()
    create_progressbar(boot, name) =
        Progress(boot, name, dt=0.1, enabled=progressbar_enabled(), color=_progressbar_color)

    trace_infos = Dict()

    boot = 10
    # 1.
    name = "Random and sparse"
    n = 10
    K = GF(2^62 + 135)
    R, x = PolynomialRing(K, [["x$i" for i in 1:n]...])
    T = 30:30
    D = 0:6
    trace_infos[(name=name, nvars=n, boot=boot, terms=T, degrees=n * D)] = []
    prog = create_progressbar(boot, "# Benchmarking.. 1/4")
    i = 1
    while i <= boot
        A = rand(R, T, D) + rand(K)
        A = A + x[1]^(total_degree(A) + 1)
        B = rand(R, T, D) + rand(K)
        B = B + x[1]^(total_degree(B) + 1)
        if length(factor(A)) > 1 || length(factor(B)) > 1
            continue
        end
        F = A * B
        trace_info = []
        A_, B_ = ExactSparseInterpolations.iterative_factor(F, trace_info=trace_info)
        push!(trace_infos[(name=name, nvars=n, boot=boot, terms=T, degrees=n * D)], trace_info)
        @assert A_ * B_ == F
        i += 1
        update!(prog, i, valuecolor=_progressbar_value_color, showvalues=[(:Iteration, i)])
    end

    # 2. 
    name = "Example by Joris"
    n = 10
    K = GF(2^62 + 135)
    R, x = PolynomialRing(K, [["x$i" for i in 1:n]...])
    R2, x2 = PolynomialRing(K, [["x$i" for i in 1:(2n)]...])
    T = 30:30
    D = 0:6
    trace_infos[(name=name, nvars=2n, boot=boot, terms=T, degrees=n * D)] = []
    prog = create_progressbar(boot, "# Benchmarking.. 2/4")
    i = 1
    while i <= boot
        A = rand(R, T, D) + rand(K)
        A = A + x[1]^(total_degree(A) + 1)
        B = rand(R, T, D) + rand(K)
        B = B + x[1]^(total_degree(B) + 1)
        A = evaluate(A, x2[1:n] .* vcat(x2[1], x2[2:n]))
        B = evaluate(B, x2[1:n] .* vcat(x2[1], x2[2:n]))
        if length(factor(A)) > 1 || length(factor(B)) > 1
            continue
        end
        F = A * B
        trace_info = []
        A_, B_ = ExactSparseInterpolations.iterative_factor(F, trace_info=trace_info)
        push!(trace_infos[(name=name, nvars=2n, boot=boot, terms=T, degrees=n * D)], trace_info)
        @assert A_ * B_ == F
        i += 1
        update!(prog, i, valuecolor=_progressbar_value_color, showvalues=[(:Iteration, i)])
    end

    # 3.
    name = "Quite dense"
    n = 10
    density = 0.3
    K = GF(2^62 + 135)
    R, x = PolynomialRing(K, [["x$i" for i in 1:n]...])
    D = 3
    trace_infos[(name=name, nvars=n, boot=boot, degrees=2D)] = []
    prog = create_progressbar(boot, "# Benchmarking.. 3/4")
    i = 1
    while i <= boot
        A = random_poly_of_density(R, D, density)
        A = A + x[1]^(total_degree(A) + 1) + rand(K)
        B = random_poly_of_density(R, D, density)
        B = B + x[1]^(total_degree(B) + 1) + rand(K)
        if length(factor(A)) > 1 || length(factor(B)) > 1
            continue
        end
        F = A * B
        trace_info = []
        A_, B_ = ExactSparseInterpolations.iterative_factor(F, trace_info=trace_info)
        push!(trace_infos[(name=name, nvars=n, boot=boot, degrees=2D)], trace_info)
        @assert A_ * B_ == F
        i += 1
        update!(prog, i, valuecolor=_progressbar_value_color, showvalues=[(:Iteration, i)])
    end

    # 4.
    name = "(Almost) embedded in a line"
    n = 10
    boot = 3
    K = GF(2^62 + 135)
    R, x = PolynomialRing(K, [["x$i" for i in 1:n]...])
    T = 10:10
    D = 100
    trace_infos[(name=name, nvars=n, boot=boot, degrees=D)] = []
    prog = create_progressbar(boot, "# Benchmarking.. 4/4")
    i = 1
    while i <= boot
        A = (x[1]^100 + sum([prod(x) * last(x)^i for i in 1:65]) + rand(K))
        B = (x[1]^99 + sum([prod(x) * last(x)^i for i in 1:55]) + rand(K))
        if length(factor(A)) > 1 || length(factor(B)) > 1
            continue
        end
        F = A * B
        trace_info = []
        A_, B_ = ExactSparseInterpolations.iterative_factor(F, trace_info=trace_info)
        push!(trace_infos[(name=name, nvars=n, boot=boot, degrees=D)], trace_info)
        @assert A_ * B_ == F
        i += 1
        update!(prog, i, valuecolor=_progressbar_value_color, showvalues=[(:Iteration, i)])
    end

    trace_infos
end

function analyze(benchmark_results)
    function aggregate(sym, results, by=mean)
        res = []
        for i in 1:length(results[1])
            push!(res, by(map(x -> getfield(x[i], sym), results)))
        end
        res
    end

    subplots = []
    for (id, results) in benchmark_results
        iterations = map(x -> x.nvars, results[1])
        len_F_mean = aggregate(:len_F, results, mean)
        sparsity_F = len_F_mean[end] / number_of_monoms_up_to_degree(id.nvars, id.degrees[end])
        len_A_mean = aggregate(:len_A, results, mean)
        len_B_mean = aggregate(:len_B, results, mean)
        len_F_std = aggregate(:len_F, results, std)
        len_A_std = aggregate(:len_A, results, std)
        len_B_std = aggregate(:len_B, results, std)
        p = plot(
            iterations,
            len_F_mean,
            yerr=len_F_std,
            linewidth=2,
            label="# terms in F",
            color=:blue,
            legend=:bottomright,
            title=string(id.name),
            titlefont=font(10),
            guidefont=font(8)
        )
        yaxis!("Number of terms")
        xaxis!("Iterations")
        # annotate!(5, 0, text("Sparsity of F: $(round(sparsity_F, digits=3))", 7, color=:blue))
        push!(subplots, p)
    end
    plot(subplots..., layout=(2, 2), size=(900, 700), margin=7mm)
end

results = benchmark()

analyze(results)
savefig((@__DIR__) * "/graphs/iterative-factor.pdf")
