
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

"""
    parent_ring_change(poly, new_ring)

Converts a polynomial to a different polynomial ring
Input
- `poly` - a polynomial to be converted
- `new_ring` - a polynomial ring such that every variable name appearing in poly appears among the generators

Output:
- a polynomial in `new_ring` “equal” to `poly`
"""
function parent_ring_change(
    poly::MPolyElem,
    new_ring::MPolyRing;
    matching = :byname,
    shift = 0,
)
    old_ring = parent(poly)
    # Construct a mapping for the variable indices.
    # Zero indicates no image of the old variable in the new ring  
    var_mapping = zeros(Int, max(nvars(old_ring), nvars(new_ring)))
    if matching === :byname
        old_symbols, new_symbols = symbols(old_ring), symbols(new_ring)
        for i in 1:length(old_symbols)
            u = old_symbols[i]
            found = findfirst(v -> (u === v), new_symbols)
            isnothing(found) && continue
            var_mapping[i] = found
        end
    elseif matching === :byindex
        var_mapping[1:(nvars(new_ring) - shift)] .= (1 + shift):nvars(new_ring)
    else
        throw(Base.ArgumentError("Unknown matching type: $matching"))
    end
    # Hoist the compatibility check out of the loop
    for i in 1:nvars(old_ring)
        if degree(poly, i) > 0 && iszero(var_mapping[i])
            throw(
                Base.ArgumentError(
                    """
                    The polynomial $poly contains a variable $(gens(old_ring)[i]) not present in the new ring.
                    New ring variables are $(gens(new_ring)))""",
                ),
            )
        end
    end
    bring = base_ring(new_ring)
    exps = Vector{Vector{Int}}(undef, length(poly))
    coefs = map(c -> bring(c), coefficients(poly))
    @inbounds for i in 1:length(poly)
        evec = exponent_vector(poly, i)
        new_exp = zeros(Int, nvars(new_ring))
        for i in 1:length(evec)
            iszero(var_mapping[i]) && continue
            new_exp[var_mapping[i]] = evec[i]
        end
        exps[i] = new_exp
    end
    return new_ring(coefs, exps)
end
