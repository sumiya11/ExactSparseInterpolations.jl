using DataFrames, StatsPlots, IndexedTables
using PrettyTables
using JLD

chunk = load((@__DIR__) * "/data_terms.jld")
data = chunk["data"]
nterms = chunk["nterms"]
nvars = chunk["nvars"]
degrees = chunk["degrees"]

formatter = (v, i, j) -> round(v, digits=3)
pretty_table(data[:, :, end], 
    row_labels=nterms, 
    header=nvars, 
    title="# terms / # variables. Runtimes (s)",
    tf=tf_markdown,
    formatters=formatter
)
println()
pretty_table(data[:, end, :], 
    row_labels=nterms, 
    header=degrees, 
    tf=tf_markdown,
    title="# terms / # degrees. Runtimes (s)",
    formatters=formatter
)

colors = palette(:tab10);
p1 = plot(nterms, data[:, end, end], 
    xlabel="# terms", ylabel="runtime, (s)",
    linewidth=3, color=colors[1],
    title="$(nvars[end]) vars, degree $(degrees[end])",
    titlefont=font(10)
);
p2 = plot(nvars, data[end, :, end],
    xlabel="# variables", #ylabel="runtime, (s)",
    linewidth=3, color=colors[2],
    title="$(nterms[end]) terms, degree $(degrees[end])",
    titlefont=font(10)
);
p3 = plot(degrees, data[end, end, :],
    xlabel="degree", # ylabel="runtime, (s)",
    linewidth=3, color=colors[3],
    title="$(nvars[end]) vars, $(nterms[end]) terms",
    titlefont=font(10)
);
gui()
display(plot(p1, p2, p3, layout=(1,3), legend=false, size=(800, 400)))

# plot!
gui()
plotlyjs()
display(surface(
    nvars, nterms, data[:, :, end],
    xlabel="# variables",
    ylabel="# terms",
    zlabel="Runtime (s)",
    title="Runtime of interpolation (degrees = $(degrees[end]))",
    size=(800, 600)
))


