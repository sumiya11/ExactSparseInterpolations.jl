
const _factorize_benchmarks = Dict{Symbol,Any}()

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
    _factorize_benchmarks[:v_transform_degrees] = NamedTuple{(:before, :after),Tuple{Int64,Int64}}[]
    _factorize_benchmarks[:v_transform_matrices] = Vector{Tuple{Int,Matrix{Int},Matrix{Int}}}()
    _factorize_benchmarks[:v_points_used] = Int[]
    ret
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
