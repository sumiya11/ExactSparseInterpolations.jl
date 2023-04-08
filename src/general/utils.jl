
const _factorize_benchmarks = Dict{Symbol, Any}()
const _data_to_vertex = Dict()
const _edge_to_data = Dict()

function dump_data()
    data1 = deepcopy(_data_to_vertex)
    data2 = deepcopy(_edge_to_data)
    empty!(_data_to_vertex)
    empty!(_edge_to_data)
    data1, data2
end
function dump_benchmarks()
    ret = deepcopy(_factorize_benchmarks)
    empty!(_factorize_benchmarks)
    _factorize_benchmarks[:t_removing_content] = Float64[]
    _factorize_benchmarks[:t_factoring_content] = Float64[]
    _factorize_benchmarks[:t_selecting_main_variable] = Float64[]
    _factorize_benchmarks[:t_finding_power_product] = Float64[]
    _factorize_benchmarks[:t_evaluating_coefficients] = Float64[]
    _factorize_benchmarks[:t_interpolating_coefficients] = Float64[]
    _factorize_benchmarks[:t_first_bivariate_factor] = Float64[]
    _factorize_benchmarks[:t_many_hensel_liftings] = Float64[]
    _factorize_benchmarks[:t_total] = Float64[]
    _factorize_benchmarks[:t_prim_and_sqfree] = Float64[]
    _factorize_benchmarks[:v_transform_degrees] = NamedTuple{(:before, :after), Tuple{Int64, Int64}}[]
    _factorize_benchmarks[:v_transform_matrices] = Vector{Tuple{Int, Matrix{Int}, Matrix{Int}}}()
    _factorize_benchmarks[:v_points_used] = Int[]
    _factorize_benchmarks[:v_tree] = Graphs.DiGraph()
    ret
end
function draw_graph(to, G, data1, data2)
    nodelabel = map(last, sort(collect(values(data1)), by=first))
    mainvars = map(last, nodelabel) 
    nodelabel = map(first, nodelabel)
    edges = collect(Graphs.edges(G))
    edgelabel = map(
        e -> Symbol(data2[Graphs.src(e), Graphs.dst(e)], ",", mainvars[Graphs.dst(e)]), 
        edges
    )
    nodesz = 3.5*length(first(nodelabel))
    nodesize = repeat([nodesz], length(nodelabel))
    plo = GraphPlot.gplot(
        G, 
        nodelabel=nodelabel,
        edgelabel=edgelabel,
        layout=GraphPlot.circular_layout,
        nodesize=nodesize,
        # edgelabeldistx=2, edgelabeldisty=2,
        edgelabelsize=1,
        nodelabelsize=1
        )
    draw(to, plo)
    plo
end
dump_benchmarks()

function addvertex_preserve!(G, data, key, additional=nothing)
    exists = haskey(data, key)
    if exists
        data[key][1]
    else
        Graphs.add_vertex!(G)
        v = Graphs.nv(G)
        data[key] = (v, additional)
        v
    end
end

function _savetree(key, value)
    global _factorize_benchmarks, _data_to_vertex, _edge_to_data
    status, main_var_idx, ha, f_deg, fi_degs = value
    G = _factorize_benchmarks[:v_tree]
    v_vertex = addvertex_preserve!(G, _data_to_vertex, ha, (f_deg, main_var_idx))
    for (h, d) in fi_degs
        u_vertex = addvertex_preserve!(G, _data_to_vertex, h, (d, main_var_idx))
        Graphs.add_edge!(G, (v_vertex, u_vertex))
        _edge_to_data[(v_vertex, u_vertex)] = status
    end
    nothing
end

@noinline function _savevalue(key, value)
    global _factorize_benchmarks
    key === :v_tree && return _savetree(key, value)
    if !haskey(_factorize_benchmarks, key)
        _factorize_benchmarks[key] = Vector{typeof(value)}(undef, 0)
    end
    push!(_factorize_benchmarks[key], value)
    nothing
end

macro saveval(flag, key, expr)
    return quote
        if $(esc(flag)) == true
            local res = $(esc(expr))
            _savevalue($key, res)
            res
        else
            $(esc(expr))
        end    
    end
end

macro savetime(flag, key, expr)
    return quote
        if $(esc(flag)) == true
            local t0 = time_ns()
            local res = $(esc(expr))
            local timepassed = (time_ns() - t0) / 1e9
            _savevalue($key, timepassed)
            res
        else
            $(esc(expr))
        end    
    end
end
