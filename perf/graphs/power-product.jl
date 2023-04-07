#=
    Power product substitutions for normalization.    
=#

# using the most recent local version
include("../../src/ExactSparseInterpolations.jl")

using Nemo, LinearAlgebra
using NamedArrays, Plots
using PlotlyJS

######################
## 4 variables

nice_color_table = ["#3D9970", "#FF4136", "#FF851B", "#6699CC"]

is_id(transform) = transform == Matrix(1I, size(transform, 1), size(transform, 2))
function run_and_extract_useful_info(f, main_var_idx, method)
    F_sub, new_vars, transform, invtransform, train_before, trailmonom, trailcoeff, attempts = ExactSparseInterpolations.find_power_product_substitution(
            f, main_var_idx, method)
    id = is_id(transform)
    td = total_degree(F_sub)
    attempts = attempts

    (is_identity=id, total_deg=td, checks=attempts)
end

methods = [:exhaustive, :lowerdiag_plus_minus_ones, 
            :random_nonzero_entry, :above_main_diagonal]
begin
    K = Nemo.GF(2^62 + 135)
    R, (x1,x2,x3,x4) = PolynomialRing(K, ["x$i" for i in 1:4])
    xs = gens(R)
    n = length(xs)

    boot = 100
    data = []
end
for (idx, i) in enumerate(2:6)
    F = (x1 + (x2^2 + x3^2 + x4^2)^i)*(x1 + (x2 + x3 + x4 + 1)^(i-1));
    main_var_idx = 1
    @info "$i: Length of F is $(length(F)), degree is $(total_degree(F))"
    push!(data, 
        (index=idx, 
        degree=total_degree(F),
        total_degrees=Dict(methods .=> map(_ -> [], methods)),
        checks=Dict(methods .=> map(_ -> [], methods))
        )
    )
    for j in 1:boot
        xsi = [rand(K) for _ in 2:n] .* (xs[2:end] .* (xs[2:end] .^ rand(Bool, n-1)))
        xsi = vcat([xs[1]], xsi)

        f = evaluate(F, xsi)

        info1 = run_and_extract_useful_info(f, main_var_idx, methods[1])
        info2 = run_and_extract_useful_info(f, main_var_idx, methods[2])
        info3 = run_and_extract_useful_info(f, main_var_idx, methods[3])
        info4 = run_and_extract_useful_info(f, main_var_idx, methods[4])

        info1.is_identity && info2.is_identity && info3.is_identity && info4.is_identity && continue
        @assert !info1.is_identity && !info2.is_identity && !info3.is_identity && !info4.is_identity

        tdfj = total_degree(f)
        push!(data[end].total_degrees[methods[1]], info1.total_deg / tdfj)
        push!(data[end].total_degrees[methods[2]], info2.total_deg / tdfj)
        push!(data[end].total_degrees[methods[3]], info3.total_deg / tdfj)
        push!(data[end].total_degrees[methods[4]], info4.total_deg / tdfj)

        push!(data[end].checks[methods[2]], info2.checks)
        push!(data[end].checks[methods[3]], info3.checks)
        push!(data[end].checks[methods[4]], info4.checks)
    end
end

begin
    totals = map(x -> length(x.total_degrees[methods[1]]), data)
    idxs = vcat(map(i -> repeat([data[i].degree], totals[i]), 1:length(data))...)
    total_degrees_1 = vcat(map(i -> data[i].total_degrees[methods[1]], 1:length(data))...)
    total_degrees_2 = vcat(map(i -> data[i].total_degrees[methods[2]], 1:length(data))...)
    total_degrees_3 = vcat(map(i -> data[i].total_degrees[methods[3]], 1:length(data))...)
    total_degrees_4 = vcat(map(i -> data[i].total_degrees[methods[4]], 1:length(data))...)
    checks_2 = vcat(map(i -> data[i].checks[methods[2]], 1:length(data))...)
    checks_3 = vcat(map(i -> data[i].checks[methods[3]], 1:length(data))...)
    checks_4 = vcat(map(i -> data[i].checks[methods[4]], 1:length(data))...)
end

function draw_1!(totals, idxs, 
        total_degrees_1, total_degrees_2, total_degrees_3, total_degrees_4)
    trace1 = box(
        y=total_degrees_1,
        x=idxs,
        name="Exhaustive search",
        marker_color=nice_color_table[1],
        boxpoints=false
    )
    trace2 = box(
        y=total_degrees_2,
        x=idxs,
        name="Increment 1 random entry",
        marker_color=nice_color_table[2],
        boxpoints=false
    )
    trace3 = box(
        y=total_degrees_3,
        x=idxs,
        name="Assign 5 random entries",
        marker_color=nice_color_table[3],
        boxpoints=false
    )
    trace4 = box(
        y=total_degrees_4,
        x=idxs,
        name="Assign above main diagonal",
        marker_color=nice_color_table[4],
        boxpoints=false
    )
    p1 = PlotlyJS.plot(
        [trace1, trace2, trace3, trace4], 
        Layout(
            title="Ratios of degree before to degree after substitution, $n variables",
            yaxis=attr(autorange=true, showgrid=true, zeroline=true,
                dtick=1, gridcolor="rgb(255, 255, 255)",
                gridwidth=1,
                zerolinecolor="rgb(255, 255, 255)",
                zerolinewidth=2),
            yaxis_title="Growth in total degree, times", 
            xaxis_title="Initial total degree",
            boxmode="group"
        )
    )    
    p1
end
function draw_2!(totals, idxs, checks_2, checks_3, checks_4)
    trace2 = box(
        y=checks_2,
        x=idxs,
        name="Increment 1 random entry",
        marker_color=nice_color_table[2],
        boxpoints=false
    )
    trace3 = box(
        y=checks_3,
        x=idxs,
        name="Assign 5 random entries",
        marker_color=nice_color_table[3],
        boxpoints=false
    )
    trace4 = box(
        y=checks_4,
        x=idxs,
        name="Assign above main diagonal",
        marker_color=nice_color_table[4],
        boxpoints=false
    )
    p2 = PlotlyJS.plot
    p2 = PlotlyJS.plot(
        [trace2, trace3, trace4], 
        Layout(
            title="Number of checks that the term is single, $n variables",
            yaxis_title="# checks", 
            xaxis_title="Initial total degree",
            boxmode="group"
        )
    )
end

plo = draw_1!(totals, idxs, 
    total_degrees_1, total_degrees_2, total_degrees_3, total_degrees_4
)
PlotlyJS.savefig(plo, "benchmark5_1.pdf")
plo = draw_2!(totals, idxs, 
    checks_2, checks_3, checks_4
)
PlotlyJS.savefig(plo, "benchmark5_2.pdf")

###################
# 5 variables
methods = [:exhaustive, :lowerdiag_plus_minus_ones, 
            :random_nonzero_entry, :above_main_diagonal]
begin
    K = Nemo.GF(2^62 + 135)
    R, (x1,x2,x3,x4,x5) = PolynomialRing(K, ["x$i" for i in 1:5])
    xs = gens(R)
    n = length(xs)

    cool_polynomials = [
        (x1 + (x2^2 + x3^2 + x4^2 + x5^2)^4)*(x1 + (x2 + x3 + x4 + x5 + 1)^3),
        (x1 + x2 + x3 + x4 + x5)^6,
        (x1 + x2 + x3 + x4 + x5) + (x1 + x2 + x3 + x4 + x5)^2 + (x1 + x2 + x3 + x4 + x5)^3 + (x1 + x2 + x3 + x4 + x5)^4,
        x1*x2*(x1 - x2 - x3)*(x2 - x3)*(x1 + x3 + x5)*(x5 - x4 - x3)*(x3 + x4)*(x1 + x3 + x4 + x5)*(x2 + x3 + x4)
    ]

    boot = 20
    data = []
end
for (idx, F) in enumerate(cool_polynomials)
    main_var_idx = 1
    @info "$idx: Length of F is $(length(F)), degree is $(total_degree(F))"
    push!(data, 
        (index=idx, 
        degree=total_degree(F),
        total_degrees=Dict(methods .=> map(_ -> [], methods)),
        checks=Dict(methods .=> map(_ -> [], methods))
        )
    )
    for j in 1:boot
        xsi = [rand(K) for _ in 2:n] .* (xs[2:end] .* (xs[2:end] .^ rand((0,0,0,0,0,0,0,0,0,0,0,0,1), n-1)))
        xsi = vcat([xs[1]], xsi)

        f = evaluate(F, xsi)

        @info "running $(methods[1]).."
        info1 = run_and_extract_useful_info(f, main_var_idx, methods[1])
        info2 = run_and_extract_useful_info(f, main_var_idx, methods[2])
        info3 = run_and_extract_useful_info(f, main_var_idx, methods[3])
        @info "running $(methods[4]).."
        info4 = run_and_extract_useful_info(f, main_var_idx, methods[4])

        info1.is_identity && info2.is_identity && info3.is_identity && info4.is_identity && continue
        @assert !info1.is_identity && !info2.is_identity && !info3.is_identity && !info4.is_identity

        tdfj = total_degree(f)
        push!(data[end].total_degrees[methods[1]], info1.total_deg / tdfj)
        push!(data[end].total_degrees[methods[2]], info2.total_deg / tdfj)
        push!(data[end].total_degrees[methods[3]], info3.total_deg / tdfj)
        push!(data[end].total_degrees[methods[4]], info4.total_deg / tdfj)

        push!(data[end].checks[methods[2]], info2.checks)
        push!(data[end].checks[methods[3]], info3.checks)
        push!(data[end].checks[methods[4]], info4.checks)
    end
end

begin
    type_of_poly = ["f * g", "(x + y)^n", "(x+y) + (x+y)^2 + ...", "(x-y)(y-z)..."]
    totals = map(x -> length(x.total_degrees[methods[1]]), data)
    idxs = vcat(map(i -> repeat([type_of_poly[i]], totals[i]), 1:length(data))...)
    total_degrees_1 = vcat(map(i -> data[i].total_degrees[methods[1]], 1:length(data))...)
    total_degrees_2 = vcat(map(i -> data[i].total_degrees[methods[2]], 1:length(data))...)
    total_degrees_3 = vcat(map(i -> data[i].total_degrees[methods[3]], 1:length(data))...)
    total_degrees_4 = vcat(map(i -> data[i].total_degrees[methods[4]], 1:length(data))...)
    checks_2 = vcat(map(i -> data[i].checks[methods[2]], 1:length(data))...)
    checks_3 = vcat(map(i -> data[i].checks[methods[3]], 1:length(data))...)
    checks_4 = vcat(map(i -> data[i].checks[methods[4]], 1:length(data))...)
end

function draw_1!(totals, idxs, 
        total_degrees_1, total_degrees_2, total_degrees_3, total_degrees_4)
    trace1 = box(
        y=total_degrees_1,
        x=idxs,
        name="Exhaustive search",
        marker_color=nice_color_table[1],
        boxpoints=false
    )
    trace2 = box(
        y=total_degrees_2,
        x=idxs,
        name="Increment 1 random entry",
        marker_color=nice_color_table[2],
        boxpoints=false
    )
    trace3 = box(
        y=total_degrees_3,
        x=idxs,
        name="Assign 5 random entries",
        marker_color=nice_color_table[3],
        boxpoints=false
    )
    trace4 = box(
        y=total_degrees_4,
        x=idxs,
        name="Assign above main diagonal",
        marker_color=nice_color_table[4],
        boxpoints=false
    )
    p1 = PlotlyJS.plot(
        [trace1, trace2, trace3, trace4], 
        Layout(
            title="Ratios of degree before to degree after substitution, $n variables",
            yaxis=attr(autorange=true, showgrid=true, zeroline=true,
                dtick=1, gridcolor="rgb(255, 255, 255)",
                gridwidth=1,
                zerolinecolor="rgb(255, 255, 255)",
                zerolinewidth=2),
            yaxis_title="Growth in total degree, times", 
            xaxis_title="Type of polynomial",
            boxmode="group"
        )
    )    
    p1
end
function draw_2!(totals, idxs, checks_2, checks_3, checks_4)
    trace2 = box(
        y=checks_2,
        x=idxs,
        name="Increment 1 random entry",
        marker_color=nice_color_table[2],
        boxpoints=false
    )
    trace3 = box(
        y=checks_3,
        x=idxs,
        name="Assign 1 random entry",
        marker_color=nice_color_table[3],
        boxpoints=false
    )
    trace4 = box(
        y=checks_4,
        x=idxs,
        name="Assign above main diagonal",
        marker_color=nice_color_table[4],
        boxpoints=false
    )
    p2 = PlotlyJS.plot
    p2 = PlotlyJS.plot(
        [trace2, trace3, trace4], 
        Layout(
            title="Number of checks that the term is single, $n variables",
            yaxis_title="# checks", 
            xaxis_title="Type of polynomial",
            boxmode="group"
        )
    )
end

plo = draw_1!(totals, idxs, 
    total_degrees_1, total_degrees_2, total_degrees_3, total_degrees_4
)
PlotlyJS.savefig(plo, "benchmark6_1.pdf")
draw_2!(totals, idxs, 
    checks_2, checks_3, checks_4
)
PlotlyJS.savefig(plo, "benchmark6_2.pdf")

#######
# 10 variables
methods = [:exhaustive, :lowerdiag_plus_minus_ones, 
            :random_nonzero_entry, :above_main_diagonal]
begin
    K = Nemo.GF(2^62 + 135)
    R, xs = PolynomialRing(K, ["x$i" for i in 1:10])
    xs = gens(R)
    n = length(xs)

    cool_polynomials = [
        (xs[1] + prod(xs[2:end]))*(xs[1] + sum(xs[2:end])^2),
        sum(xs) + sum(xs)^2 + sum(xs)^3
    ]

    boot = 100
    data = []
end
for (idx, F) in enumerate(cool_polynomials)
    main_var_idx = 1
    @info "$idx: Length of F is $(length(F)), degree is $(total_degree(F))"
    push!(data, 
        (index=idx, 
        degree=total_degree(F),
        total_degrees=Dict(methods .=> map(_ -> [], methods)),
        checks=Dict(methods .=> map(_ -> [], methods))
        )
    )
    for j in 1:boot
        xsi = [rand(K) for _ in 2:n] .* (xs[2:end] .* (xs[2:end] .^ rand((0,0,0,0,0,0,0,0,0,0,0,0,1), n-1)))
        xsi = vcat([xs[1]], xsi)

        f = evaluate(F, xsi)

        info2 = run_and_extract_useful_info(f, main_var_idx, methods[2])
        info3 = run_and_extract_useful_info(f, main_var_idx, methods[3])
        @info "running $(methods[4]).."
        info4 = run_and_extract_useful_info(f, main_var_idx, methods[4])

        info2.is_identity && info3.is_identity && info4.is_identity && continue
        @assert !info2.is_identity && !info3.is_identity && !info4.is_identity

        tdfj = total_degree(f)
        push!(data[end].total_degrees[methods[2]], info2.total_deg / tdfj)
        push!(data[end].total_degrees[methods[3]], info3.total_deg / tdfj)
        push!(data[end].total_degrees[methods[4]], info4.total_deg / tdfj)

        push!(data[end].checks[methods[2]], info2.checks)
        push!(data[end].checks[methods[3]], info3.checks)
        push!(data[end].checks[methods[4]], info4.checks)
    end
end

begin
    type_of_poly = ["f * g", "(x+y) + (x+y)^2 + ..."]
    totals = map(x -> length(x.total_degrees[methods[2]]), data)
    idxs = vcat(map(i -> repeat([type_of_poly[i]], totals[i]), 1:length(data))...)
    total_degrees_2 = vcat(map(i -> data[i].total_degrees[methods[2]], 1:length(data))...)
    total_degrees_3 = vcat(map(i -> data[i].total_degrees[methods[3]], 1:length(data))...)
    total_degrees_4 = vcat(map(i -> data[i].total_degrees[methods[4]], 1:length(data))...)
    checks_2 = vcat(map(i -> data[i].checks[methods[2]], 1:length(data))...)
    checks_3 = vcat(map(i -> data[i].checks[methods[3]], 1:length(data))...)
    checks_4 = vcat(map(i -> data[i].checks[methods[4]], 1:length(data))...)
end

function draw_1!(totals, idxs, 
        total_degrees_2, total_degrees_3, total_degrees_4)
    trace2 = box(
        y=total_degrees_2,
        x=idxs,
        name="Increment 1 random entry",
        marker_color=nice_color_table[2],
        boxpoints=false
    )
    trace3 = box(
        y=total_degrees_3,
        x=idxs,
        name="Assign 5 random entries",
        marker_color=nice_color_table[3],
        boxpoints=false
    )
    trace4 = box(
        y=total_degrees_4,
        x=idxs,
        name="Assign above main diagonal",
        marker_color=nice_color_table[4],
        boxpoints=false
    )
    p1 = PlotlyJS.plot(
        [trace2, trace3, trace4], 
        Layout(
            title="Ratios of degree before to degree after substitution, $n variables",
            yaxis=attr(autorange=true, showgrid=true, zeroline=true,
                dtick=3, gridcolor="rgb(255, 255, 255)",
                gridwidth=1,
                zerolinecolor="rgb(255, 255, 255)",
                zerolinewidth=2),
            yaxis_title="Growth in total degree, times", 
            xaxis_title="Type of polynomial",
            boxmode="group"
        )
    )    
    p1
end
function draw_2!(totals, idxs, checks_2, checks_3, checks_4)
    trace2 = box(
        y=checks_2,
        x=idxs,
        name="Increment 1 random entry",
        marker_color=nice_color_table[2],
        boxpoints=false
    )
    trace3 = box(
        y=checks_3,
        x=idxs,
        name="Assign 1 random entry",
        marker_color=nice_color_table[3],
        boxpoints=false
    )
    trace4 = box(
        y=checks_4,
        x=idxs,
        name="Assign above main diagonal",
        marker_color=nice_color_table[4],
        boxpoints=false
    )
    p2 = PlotlyJS.plot
    p2 = PlotlyJS.plot(
        [trace2, trace3, trace4], 
        Layout(
            title="Number of checks that the term is single, $n variables",
            yaxis_title="# checks", 
            xaxis_title="Type of polynomial",
            boxmode="group"
        )
    )
end

plo = draw_1!(totals, idxs, 
    total_degrees_2, total_degrees_3, total_degrees_4
)
PlotlyJS.savefig(plo, "benchmark7_1.pdf")
draw_2!(totals, idxs, 
    checks_2, checks_3, checks_4
)
# PlotlyJS.savefig(plo, "benchmark5_2.pdf")



# speedup = NamedArray(round.(times1./times2,sigdigits=3))
# setnames!(speedup, string.(packabletypes), 1)
# setnames!(speedup, string.(packabletypes), 2)
# println("\nspeedup:")
# speedup

