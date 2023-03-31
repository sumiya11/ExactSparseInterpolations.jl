# Factorization!

##########################
##########################

# Returns the index of the main variable
# (or 0, signalling that no such variable exists)
#
# Does n bivariate factorizations at worst
function select_main_variable(F)
    ring = parent(F)
    K, xs = base_ring(ring), gens(ring)
    biv_ring, (x, u) = PolynomialRing(K, ["x","u"])
    univ_ring, t = PolynomialRing(K, "t")
    # traverse variables in the order of increasing degree
    var2degree = [(i, degree(F, xs[i])) for i in 1:length(xs)]
    sort!(var2degree, by=x->x[2])
    for (i, deg) in var2degree
        deg < 2 && continue
        p = vcat(zeros(ring, i - 1), [xs[i]], zeros(ring, length(xs) - i))
        lc = leading_coefficient_in(F, xs[i])
        # F(0, .., 0, xi, 0, .., 0) is squarefree
        iszero(evaluate(lc, p)) && continue
        !issquarefree(to_univariate(univ_ring, evaluate(F, p))) && continue
        # F(u, .., u, xi, u, .., u) is decomposes into > 1 factors
        up = vcat(repeat([u], i - 1), [x], repeat([u], length(xs) - i))
        F_u = evaluate(F, up)
        Fi, fi, Si = revealing_bivariate_factorization_ff(F_u)
        length(Fi) < 2 && continue
        # 
        return i
    end
    0
end

# Checks that f has a single term independent of the main variable
# that has the smallest total degree
function check_invariant(f, main_var_idx)
    monoms = filter(m -> degree(m, main_var_idx) == 0, collect(monomials(f)))
    min_deg = minimum(total_degree, monoms)
    cnt = count(m -> total_degree(m) == min_deg, monoms)
    trailmonom = monoms[findfirst(m -> degree(m, main_var_idx) == 0 && total_degree(m) == min_deg, monoms)]
    cnt == 1, trailmonom
end

function transform_poly(f, transform)
    ring = parent(f)
    monoms = collect(monomials(f))
    newmonoms = map(m -> transform * exponent_vector(m, 1), monoms)
    ring(collect(coefficients(f)), newmonoms)
end

function next_monom_transform(transform, main_var_idx, monomtransform)
    if monomtransform === :lowerdiag_plus_minus_ones
        generate_next_transform_1(transform, main_var_idx)
    elseif monomtransform === :product_of_unimodular
        generate_next_transform_2(transform, main_var_idx)
    else
        generate_next_transform_1(transform, main_var_idx)
    end
end

# Slighly change the values under the main diagonal
function generate_next_transform_1(transform, main_var_idx)
    n = size(transform, 1)
    i, j = 1, 1
    while i == j || i == main_var_idx || j == main_var_idx
        i = rand(1:n)
        j = rand(1:max(1, i))
    end
    # I just really don't like minus ones :)
    transform[i, j] += 1
    transform
end

# Multiply by a unimodular matrix
function generate_next_transform_2(transform, main_var_idx)
    n = size(transform, 1)
    i, j = 1, 1
    while i == j || i == main_var_idx || j == main_var_idx
        i = rand(1:n)
        j = rand(i:n)
    end
    U = Matrix(1I, n, n)
    if rand(Bool)
        j, i = i, j
    end
    U[i, j] = 1
    transform = transform * U
    transform
end

# Multiply by
# [ 1  +1  0    0 ]
# [ 0   1  +1   0 ]
# [ 0   0   1   +1]
# [ 0   0   0    1 ] 
function generate_next_transform_3(transform, main_var_idx)
    n = size(transform, 1)
    U = Matrix(1I, n, n)
    for i in 1:n-1
        if i == main_var_idx
            continue
        end
        U[i, i+1] = 1 * rand(Bool)
    end
    U
end

# Returns the 
# (f_subs, subs, transform, invtransform, trailmonom)
# . f_subs - transformed polynomial
# . subs - variable exponents for the substitution
# . transform - the matrix of the transformation in the exponent space 
# . trailmonom
function find_power_product_substitution(f, main_var_idx, monomtransform)
    # There must be a single term independent of the main variable 
    # and having with the lowest total degree..
    ring = parent(f)
    xs = gens(ring)
    n = length(xs)
    # ..find the transform on the degrees 
    transform = Matrix(1I, n, n)
    new_vars = xs
    f_subs, trailmonom = f, f
    k, attempts = 0, 100
    while k < attempts
        k += 1
        # ..obtain the transformed polynomial
        new_vars = map(var -> transform_poly(var, transform), xs)
        f_subs = evaluate(f, new_vars) 
        success, trailmonom = check_invariant(f_subs, main_var_idx)
        success && break
        # ..scramble the transform a bit
        transform = next_monom_transform(transform, main_var_idx, monomtransform)
    end
    # he-he
    invtransform = round.(Int, inv(transform))
    f_subs, new_vars, transform, invtransform, trailmonom, attempts
end

# Evaluate the d*D coefficients of the polynomial F_u at the point c.
function evaluate_coefficients(F_u, var_main, var_u, c)
    point = vcat([var_main, var_u], c)
    evaluate(F_u, point)
end

# Evaluate the d*D coefficients of the polynomial F_u 
# at the sequence of points ci.
function evaluate_coefficients(F_u, var_main, var_u, ci::Vector{Vector{T}}) where {T}
    map(c -> evaluate_coefficients(F_u, var_main, var_u, c), ci)
end

function fill_coeff_evaluations!(
        Fcoeffs::Vector{Vector{Vector{T}}}, Fi, 
        k::Integer, npoints::Integer) where {T}
    K = base_ring(parent(Fi[1]))
    if !isassigned(Fcoeffs, 1)
        for i in 1:length(Fi)
            Fcoeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(Fi[i]))
            for j in 1:length(Fcoeffs[i])
                Fcoeffs[i][j] = Vector{elem_type(K)}(undef, npoints)
            end
        end
    end
    @inbounds for i in 1:length(Fi)
        for j in 1:length(Fcoeffs[i])
            Fcoeffs[i][j][k] = Nemo.coeff(Fi[i], j)
        end
    end
    nothing
end

function normalize_factorization(Fi::Vector{T}, trail) where {T}
    ai = map(trailing_coefficient, Fi)
    Fi = map((f, c) -> trail*divexact(f, c), Fi, ai)
    Fi
end

# Returns a return code and potential factors of F.
# Assumes there are <= T terms to interpolate at least in some of the factors.
function _find_some_factors(
            F, T::Integer,
            bench::Bool,
            monomtransform::Symbol
        )
    @assert T > 0
    ring = parent(F)
    K, xs = base_ring(ring), gens(ring)
    n = length(xs)
    @assert characteristic(K) > 0

    # Select the main variable
    @savetime bench :t_select_main_variable main_var_idx = select_main_variable(F)
    @saveval bench :v_main_var (main_var_idx=main_var_idx, degrees=map(i -> degree(F, i), 1:n))
    # @info "Main variable" main_var_idx F T
    # No good main variable has been found
    if main_var_idx == 0
        return false, F
    end
    
    # Find the variable substitution that gives a normalization
    F_sub, new_vars, transform, invtransform, trailmonom, attempts = @savetime bench :t_find_power_product find_power_product_substitution(F, main_var_idx, monomtransform)
    @saveval bench :v_transform_matrices (attempts, transform, invtransform)
    @saveval bench :v_transform_degrees (total_degree(F), total_degree(F_sub))
    # just to be sure
    trailmonom = trailmonom

    # @info "Transforms"  transform invtransform
    d_main, D_total = degree(F, main_var_idx), total_degree(F)
    # @info "Degrees before subst.:" d_main D_total
    d_main, D_total = degree(F_sub, main_var_idx), total_degree(F_sub)
    # @info "Degrees after subst.:" d_main D_total

    # Initialize the interpolation routine
    npoints = 2T
    ring_coeff, _ = PolynomialRing(K, ["y$i" for i in 2:n])
    ring_u, (t_u, u_u) = PolynomialRing(K, ["t", "u"])
    ring_ext, (t_ext, u_ext, x_ext...) = PolynomialRing(ring_u, vcat(["t", "u"], ["y$i" for i in 2:n]))
    partial_degrees = [degree(F_sub, i) + degree(trailmonom, i) for i in 1:n if i != main_var_idx]
    interpolator = FasterMultivariateBenOrTiwari(ring_coeff, T, partial_degrees)
    # cs = [c2...cn] is an element of K^(n-1), which will be the starting
    # point of the geometric sequence as requested by the interpolator
    cs = startingpoint(interpolator)
    
    # Reduce F to a bivariate polynomial in the main variable and u:
    # evaluate F at (x1, x2 u, x3 u, ..., xn u)
    point_u = vcat(
        u_ext .* x_ext[1:main_var_idx-1], 
        t_ext,
        u_ext .* x_ext[main_var_idx:end]
    )
    F_sub_u = evaluate(F_sub, point_u)
    # evaluate F at (t, c1^0 u, c2^0 u, ..., cn^0 u)
    # (the first point in the sequence will be c^0).
    F_sub_u_c = evaluate_coefficients(F_sub_u, t_u, u_u, cs .^ 0)
    trail = evaluate(trailmonom, vcat([one(K)], cs .^ 0))
    
    # @info "Evaluated polys" F_sub_u F_sub_u_c trail typeof(F_sub_u_c)
    # Factor bivariate polynomial in K[t, u].
    #
    # Reveals the bivariate factors Fi, together with
    # univariate factors f1,...,fk and disjoint sets Si, so that
    #   F_c = prod_i Fi
    #   F_c(t, u) = f1(t)*...*fk(t)  mod(u) 
    #   Fi(t, u) = prod_i f_i(t)  mod(u), i in Si
    @savetime bench :t_first_bivariate_factor Fi, fi, Si = revealing_bivariate_factorization_ff(F_sub_u_c)
    # @info "Revealing factorization" Fi fi Si typeof(Fi) typeof(fi)
    Fi = normalize_factorization(Fi, trail)
    
    # The storages to collect the coefficients of factors 
    # for different evaluation points.
    # Initialize storages and fill the values for the first evaluation point.
    Fcoeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(Fi))
    fill_coeff_evaluations!(Fcoeffs, Fi, 1, npoints)
    
    # Evaluate another 2T - 1 times
    cpoints = map(i -> cs .^ i, 0:npoints-1)
    @savetime bench :t_evaluating_coefficients F_sub_u_ci = evaluate_coefficients(F_sub_u, t_u, u_u, cpoints)
    @saveval bench :v_points_used length(cpoints)

    @savetime bench :t_many_hensel_liftings for j in 2:npoints
        F_sub_u_c = F_sub_u_ci[j]
        # Here, all of these factorizations should benefit from the fact
        # that we already know correct in univariate factorization
        # and the corresponding recombination pattern.
        # We just need to hensel lift to F_c(t, u) mod <u^D>
        Fi = bivariate_factorization_ff_with_known_univariate(F_sub_u_c, fi, Si)
        trail = evaluate(trailmonom, vcat([one(K)], cpoints[j]))
        Fi = normalize_factorization(Fi, trail)
        # Collect the coefficients
        fill_coeff_evaluations!(Fcoeffs, Fi, j, npoints)
    end
    # Interpolate !
    coeffs_interpolated = Vector{Vector{elem_type(ring_coeff)}}(undef, length(Fcoeffs))
    @savetime bench :t_interpolation for i in 1:length(Fcoeffs)
        coeffs_interpolated[i] = Vector{elem_type(ring_coeff)}(undef, length(Fcoeffs[i]))
        for j in 1:length(Fcoeffs[i])
            coeffs_interpolated[i][j] = interpolate!(interpolator, cpoints, Fcoeffs[i][j])
        end
    end
    
    # @info "Interpolated coeffs" coeffs_interpolated
    # Select only factors which have less than a half of the coefficients present.
    # This means that 
    #=
        number of interpolated terms < T/2
    or, on average:
        number of interpolated terms < T/2 dD
    =#
    # Check 1.5 instead of 2
    good_factors = filter(
        i -> all(c -> length(c) <= div(T, 2), coeffs_interpolated[i]), 
        1:length(Fcoeffs)
    )
    if isempty(good_factors)
        return false, Vector{typeof(F)}(undef, 0)
    end
    # Collect the interpolated coefficients into a polynomial.
    Fmonoms = map(collect ∘ monomials, Fi)
    Finterpolated = map(
        (ms, cfs) -> 
            sum(evaluate(tt, [xs[main_var_idx], one(ring)])*evaluate(cc, xs[1:n .!= main_var_idx]) 
            for (tt, cc) in zip(ms, cfs)), 
        Fmonoms[good_factors],
        coeffs_interpolated[good_factors]
    )
    # @info "after unify" Finterpolated
    # perform backwards substitution
    Finterpolated = map(f -> transform_poly(f, invtransform), Finterpolated)
    # @info "after inv subst" Finterpolated
    # divide out the content in the main variable
    Finterpolated = map(f -> divexact(f, content_in(f, xs[main_var_idx])), Finterpolated)
    # @info "after remove content" Finterpolated
    # normalize by the constant leading coefficient
    Finterpolated = map(f -> divexact(f, leading_coefficient(f)), Finterpolated)
    
    true, Finterpolated
end

function _factorize_recursive_prim_squarefree(
        F, bench, monomtransform; 
        ubT=2length(F))  # 1. this is a heurusitic!
    @assert ubT > 0
    T = 2
    success, Pi = false, [F]
    @label Start
    while T <= ubT
        success, Pi = _find_some_factors(F, T, bench, monomtransform)
        T = 2*T
        success && break
    end
    !success && return [F]
    P = prod(Pi)
    flag, Q = divides(F, P)  # 2. this is a heurusitic!
    flag && isunit(Q) && return Pi
    !flag && @goto Start
    Qi = _factorize_recursive_prim_squarefree(
        Q, bench, monomtransform, ubT=ubT)
    return vcat(Pi, Qi)
end

function _factorize(
        F, 
        strategy, normalization, mainvar, bench, 
        skipcontent, monomtransform, interpolator)
    # We assume F is squarefree --
    # there is no squarefree factorization implemented :^(
    
    # Get the primitive part of F.
    # Don't forget to recursively factorize the content
    if skipcontent
        F_prim = F
    else
        @savetime bench :t_removing_content F_prim = my_primpart(F)
    end

    # Factor the primitive part.
    if strategy === :recursive
        Fi = _factorize_recursive_prim_squarefree(F_prim, bench, monomtransform)
    else
        Fi = _factorize_recursive_prim_squarefree(F_prim, bench, monomtransform)
    end
    
    # Factor the content.
    if !skipcontent
        @savetime bench :t_factoring_content begin
            cont = divexact(F, F_prim)
            if !isone(cont)
                Fj = _factorize(
                    cont, strategy, normalization, mainvar, 
                    bench, skipcontent, monomtransform, interpolator)
                append!(Fi, Fj)
            end
        end
    end

    Fi
end

##########################
##########################

_supported_keywords = [
    :strategy, :normalization, :mainvar, :benchmark,
    :skipcontent, :monomtransform, :interpolator
]

# Top level factorization.
# Supported keyword arguments are:
# - `strategy`, either `:recursive` or `:revealing`,
# - `normalization`, either `:shift` or `:powerproduct`,
# - `monomtransform`, either of `:lowerdiag_plus_minus_ones`, `:product_of_unimodular`
#   `:lowerdiag_plus_minus_ones` will generate matrices with +-1 below the main diagonal 
#   with the increasing density.
#   `:product_of_unimodular` will generate matrices which are products of simple 
#   unimodular upper and lower triangular matrices. 
# - `mainvar`, either `:random`, or `:smalldegree`, or `:bigdegree`,
# - `benchmark`, either `true` or `false`,
# - `skipcontent`, either `true` or `false`:
#   use `true` if you know that the content is 1.
# - `interpolator`, possible options are
#   `:automatically`, `:kron_benortiwari`, :`prime_benortiwari`
function top_level_factorize(F; kws...)
    strategy = get(kws, :strategy, :recursive) # default strategy is recursive
    normalization = get(kws, :normalization, :powerproduct) # default normalization is powerproduct
    mainvar = get(kws, :mainvar, :smalldegree) # default main variable is smalldegree
    benchmark = get(kws, :benchmark, false) # default benchmark is false
    skipcontent = get(kws, :skipcontent, false) # default skipcontent is false
    monomtransform = get(kws, :monomtransform, :lowerdiag_plus_minus_ones) # default monomtransform is lowerdiag_plus_minus_ones
    interpolator = get(kws, :interpolator, :automatically) # default interpolator is automatically
    if any(key -> !(key ∈ _supported_keywords), keys(kws))
        @warn "Some keyword arguments are not not supported. Supported ones are:\n$_supported_keywords"
    end
    @savetime benchmark :t_total Fi = _factorize(
        F, strategy, normalization, mainvar, benchmark,
        skipcontent, monomtransform, interpolator
    )
    Fi
end
