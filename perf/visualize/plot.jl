using DataFrames, StatsPlots, IndexedTables

# load from the file
data = load((@__DIR__) * "/data.jld")
data = data["data"]

# plot!
plotlyjs()
surface(
    nterms, nvars, data[:, :, end],
    xlabel="# terms",
    ylabel="# variables",
    zlabel="Runtime (s)",
    title="Runtime of interpolation (degrees = 32)",
    size=(800, 600)
)
