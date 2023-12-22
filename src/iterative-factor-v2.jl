
function iterative_factor(F; trace_info=[])
    @assert length(Nemo.factor(F)) == 2  # just to be sure!
    ring = parent(F)
    n = nvars(ring)
    field = base_ring(ring)
    c, F2, bivariate_factors = nothing, nothing, nothing
    # the projection may split into 3+ factors, want to avoid that
    while true
        c = distinct_points(field, n)
        # c = map(field, [i + 6 for i in 1:n])
        F2 = project(F, 2, c)
        bivariate_factors = collect(Nemo.factor(F2))
        length(bivariate_factors) == 2 && break
    end
    @info "" F2 bivariate_factors
    A2, B2 = map(first, collect(Nemo.factor(F2)))
    A2 = leading_coefficient(F2) * A2
    B2 = leading_coefficient(F2) * B2
    @info """
    Projection point: $c
    Bivariate factors: 
    \t$A2
    \t$B2"""
    An, Bn = _iterative_factor(F, A2, B2, c, trace_info=trace_info)
    An, Bn
end

function _iterative_factor(F, A2, B2, c; trace_info=[])
    ring = parent(F)
    field = base_ring(ring)
    n = nvars(ring)
    # Base case:
    if n == 2
        push!(trace_info, (nvars=n, len_A=length(A2), len_B=length(B2), len_F=length(F)))
        return A2, B2
    end
    # Recursive call to obtain the factorization in n - 1 variables.
    k = n - 1
    Fk = project(F, k, c)
    Ak, Bk = _iterative_factor(Fk, A2, B2, c, trace_info=trace_info)
    if length(Ak) > length(Bk)
        # Interpolate the smallest factor
        Ak, Bk = Bk, Ak
    end
    @info "Lifting from $(n - 1) to $n variables"
    @info """
    Fk = $Fk
    Ak = $Ak
    Bk = $Bk"""
    # Evaluating F and Ak at a geometric progression.
    # (assuming T points are enough)
    T = 2length(F)
    point = distinct_points(field, k)
    # point = map(field, Primes.primes(2, 100)[1:k])
    @info """
    Generator: $point
    Evaluation points: $T"""
    # Densification in t for Ak(x1 t, ..., xk t, ck+1, ..., cn) and analogous for Bk
    Ak_t = apply_regularizing_weight(Ak, ones(Int, k))
    Bk_t = apply_regularizing_weight(Bk, ones(Int, k))
    Ak_t_coeffs_eval = simultaneous_multivariate_evaluate(ring, Ak_t, point, T)
    Bk_t_coeffs_eval = simultaneous_multivariate_evaluate(ring, Bk_t, point, T)
    # First, densification in u for F(x1 t, ..., xk t, u)
    F_u = apply_regularizing_weight(F, vcat(zeros(Int, k), 1))
    # Second, densification in t for each u-densified chunk of F
    F_tu = map(f -> apply_regularizing_weight(f, vcat(ones(Int, k), 0)), F_u)
    F_tu_coeffs_eval = map(
        f_coeffs ->
            simultaneous_multivariate_evaluate(ring, f_coeffs, vcat(point, one(field)), T),
        F_tu
    )
    @info """
    Densification in t:
    \tAk_t = $Ak_t
    \tBk_t = $Bk_t
    Densification in t and u (variable $(gen(ring, n)) maps to u):
    \tF_u = $F_u
    \tF_tu = $F_tu"""
    # Hensel lift Ak and Ak to An and An.
    # Note the modulo used in hensel lifting.
    ring_bivariate, (t, u) = PolynomialRing(field, ["t", "u"])
    degree_in_u = length(F_tu)
    lifting_bound = degree_in_u + 1
    lifting_modulo = u - c[n]
    A_tu_lifted = Vector{typeof(Ak)}(undef, T)
    B_tu_lifted = Vector{typeof(Ak)}(undef, T)
    @debug """
    Hensel lifting bound: O(($lifting_modulo)^$lifting_bound)
    Hensel lifting modulo: $lifting_modulo"""
    for i in 1:T
        # Here, I dance with a tambourine to instantiate the evaluations of Ak,
        # Bk, and F as the polynomials acceptable by the Hensel lifting routine
        Ak_t_i =
            ring_bivariate(Ak_t_coeffs_eval[i], [[j - 1, 0] for j in 1:length(Ak_t_coeffs_eval[i])])
        Bk_t_i =
            ring_bivariate(Bk_t_coeffs_eval[i], [[j - 1, 0] for j in 1:length(Bk_t_coeffs_eval[i])])
        F_i_coeffs = reduce(vcat, [F_tu_coeffs_eval[j][i] for j in 1:length(F_tu_coeffs_eval)])
        F_i_exps = [
            [j - 1, k - 1] for k in 1:length(F_tu_coeffs_eval) for
            j in 1:length(F_tu_coeffs_eval[k][i])
        ]
        F_tu_i = ring_bivariate(F_i_coeffs, F_i_exps)
        # NOTE: we pass monic polynomials to Hensel lifting, and multiply by the
        # corresponding (numeric) leading coefficients afterwards
        Ak_t_i_ = divexact(Ak_t_i, leading_coefficient_in(Ak_t_i, t))
        Bk_t_i_ = divexact(Bk_t_i, leading_coefficient_in(Bk_t_i, t))
        AB_tu_i = hensel_lift(F_tu_i, [Ak_t_i_, Bk_t_i_], lifting_bound, lifting_modulo)
        A_tu_lifted[i] = leading_coefficient_in(Ak_t_i, t) * AB_tu_i[1]
        B_tu_lifted[i] = leading_coefficient_in(Bk_t_i, t) * AB_tu_i[2]
    end
    @info """
    Lifted evaluated factors:
    A: $A_tu_lifted
    B: $B_tu_lifted"""
    # Interpolate the coefficient C in front of each term C t^e1 u^e2 in the
    # lifted factors. A call to interpolation is just a single solve of
    # transposed Vandermonde matrix
    ring_univariate, _ = PolynomialRing(field, "t")
    original_vars = gens(ring)
    last_var = original_vars[n]
    An, Bn = zero(ring), zero(ring)
    for i in 1:length(A_tu_lifted[1])
        kth_monom = monomial(A_tu_lifted[1], i)
        kth_term_evaluated = map(poly_eval -> coeff(poly_eval, i), A_tu_lifted)
        t_deg, u_deg = degree(kth_monom, t), degree(kth_monom, u)
        kth_term_exps = filter(monom -> total_degree(monom) == t_deg, collect(monomials(Ak)))
        kth_term_exps_evaluated = map(monom -> evaluate(monom, point), kth_term_exps)
        kth_term_evaluated = kth_term_evaluated[1:length(kth_term_exps_evaluated)]
        kth_term_coeffs = solve_transposed_vandermonde(
            ring_univariate,
            kth_term_exps_evaluated,
            kth_term_evaluated
        )
        interpolated_sum_of_terms = sum(kth_term_coeffs .* kth_term_exps)
        interpolated_sum_of_terms = evaluate(interpolated_sum_of_terms, original_vars[1:k])
        interpolated_sum_of_terms = interpolated_sum_of_terms * last_var^(u_deg)
        An += interpolated_sum_of_terms
    end
    @info "" An
    Bn = divexact(F, An)
    @info """
    The number of terms in the factors:
    An, Bn, Fn: $(length(An)), $(length(Bn)), $(length(F))
    Points used: $T"""
    push!(trace_info, (nvars=n, len_A=length(An), len_B=length(Bn), len_F=length(F)))
    An, Bn
end

####################
####################

# For example, you can run the following commands in Julia.

# # First, install / update the packages:
# using Pkg
# Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
# Pkg.add("Nemo")
#
# # Then, you should be able to do:
# using Nemo #, ExactSparseInterpolations

# R, (x1, x2, x3) = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:3]...])

# f = (x1 + 1) * (x1 + x1 * x3 + 100);
# A, B = iterative_factor(f);

# A, B

# @assert A * B == f
