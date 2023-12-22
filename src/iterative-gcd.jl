
function iterative_gcd(P, Q; trace_info=[])
    ring = parent(P)
    n = nvars(ring)
    field = base_ring(ring)
    # dilation = distinct_points(field, n)
    # P = evaluate(P, gens(ring) .* dilation)
    # Q = evaluate(Q, gens(ring) .* dilation)
    c = distinct_points(field, n)
    # c = map(field, [i + 6 for i in 1:n])
    P2, Q2 = project(P, 2, c), project(Q, 2, c)
    G2 = gcd(P2, Q2)
    @debug """
    Projection point: $c
    Bivariate polynomials: 
    \tP2 = $P2
    \tQ2 = $Q2
    \tgcd(P2, Q2) = $G2"""
    Gn = _iterative_gcd(P, Q, G2, c, trace_info=trace_info)
    # Gn = evaluate(Gn, gens(ring) .* inv.(dilation))
    Gn = divexact(Gn, leading_coefficient(Gn))
    Gn
end

function _iterative_gcd(P, Q, G2, c; trace_info=[])
    ring = parent(P)
    field = base_ring(ring)
    n = nvars(ring)
    # Base case:
    if n == 2
        # push!(trace_info, (nvars=n, len_A=length(A2), len_B=length(B2), len_F=length(F)))
        return G2
    end
    # Recursive call to obtain the gcd in n - 1 variables.
    k = n - 1
    Pk = project(P, k, c)
    Qk = project(Q, k, c)
    Gk = _iterative_gcd(Pk, Qk, G2, c, trace_info=trace_info)
    @info "Lifting from $(n - 1) to $n variables"
    @info """
    Pk = $Pk
    Qk = $Qk
    Gk = $Gk"""
    T = length(Gk) * (max(total_degree(P), total_degree(Q)) + 1)
    # point = distinct_points(field, k)
    point = map(field, Primes.primes(2, 100)[1:k])
    @debug "Evaluation points: $T"
    Gk_t = apply_regularizing_weight(Gk, ones(Int, k))
    Gk_t_coeffs_eval = simultaneous_multivariate_evaluate(ring, Gk_t, point, T)
    # First, densification in u for F(x1 t, ..., xk t, u)
    P_u = apply_regularizing_weight(P, vcat(zeros(Int, k), 1))
    Q_u = apply_regularizing_weight(Q, vcat(zeros(Int, k), 1))
    # Second, densification in t for each u-densified chunk of F
    P_tu = map(f -> apply_regularizing_weight(f, vcat(ones(Int, k), 0)), P_u)
    P_tu_coeffs_eval =
        map(cfs -> simultaneous_multivariate_evaluate(ring, cfs, vcat(point, one(field)), T), P_tu)
    Q_tu = map(f -> apply_regularizing_weight(f, vcat(ones(Int, k), 0)), Q_u)
    Q_tu_coeffs_eval =
        map(cfs -> simultaneous_multivariate_evaluate(ring, cfs, vcat(point, one(field)), T), Q_tu)
    @debug """
    Densification in t:
    \tGk_t = $Gk_t
    Densification in t and u (variable $(gens(ring)[n]) maps to u):
    \tP_tu = $P_tu
    \tQ_tu = $Q_tu"""
    # Hensel lift Ak and Ak to An and An.
    # Note the modulo used in hensel lifting.
    ring_bivariate, (t, u) = PolynomialRing(field, ["t", "u"])
    degree_in_u = min(length(P_tu), length(Q_tu))
    G_tu_lifted = Vector{typeof(Gk)}(undef, T)
    for i in 1:T
        # Here, I dance with a tambourine to instantiate the evaluations of Ak,
        # Bk, and F as the polynomials acceptable by the Hensel lifting routine
        Gk_t_i =
            ring_bivariate(Gk_t_coeffs_eval[i], [[j - 1, 0] for j in 1:length(Gk_t_coeffs_eval[i])])
        P_i_coeffs = reduce(vcat, [P_tu_coeffs_eval[j][i] for j in 1:length(P_tu_coeffs_eval)])
        P_i_exps = [
            [j - 1, k - 1] for k in 1:length(P_tu_coeffs_eval) for
            j in 1:length(P_tu_coeffs_eval[k][i])
        ]
        P_tu_i = ring_bivariate(P_i_coeffs, P_i_exps)
        Q_i_coeffs = reduce(vcat, [Q_tu_coeffs_eval[j][i] for j in 1:length(Q_tu_coeffs_eval)])
        Q_i_exps = [
            [j - 1, k - 1] for k in 1:length(Q_tu_coeffs_eval) for
            j in 1:length(Q_tu_coeffs_eval[k][i])
        ]
        Q_tu_i = ring_bivariate(Q_i_coeffs, Q_i_exps)
        # NOTE: we pass monic polynomials to Hensel lifting, and multiply by the
        # corresponding (numeric) leading coefficients afterwards
        # NOTE: Computing the full bivariate gcd here!
        G_tu_i = Nemo.gcd(P_tu_i, Q_tu_i)
        @debug "" P_tu_i Q_tu_i G_tu_i
        # G_tu_i = G_tu_i * leading_coefficient(P_tu_i) * leading_coefficient(Q_tu_i)
        G_tu_i = G_tu_i * 1
        # G_tu_i = divexact(G_tu_i, trailing_coefficient(G_tu_i))
        G_tu_lifted[i] = G_tu_i * leading_coefficient_in(Gk_t_i, t)
    end
    @info """
    Lifted evaluated gcds:
    $G_tu_lifted"""
    # Interpolate the coefficient C in front of each term C t^e1 u^e2 in the
    # lifted factors. A call to interpolation is just a single solve of
    # transposed Vandermonde matrix
    ring_univariate, _ = PolynomialRing(field, "t")
    original_vars = gens(ring)
    last_var = original_vars[n]
    Gn = zero(ring)
    for i in 1:length(G_tu_lifted[1])
        kth_monom = monomial(G_tu_lifted[1], i)
        kth_term_evaluated = map(poly_eval -> coeff(poly_eval, i), G_tu_lifted)
        @warn "" i kth_monom kth_term_evaluated
        t_deg, u_deg = degree(kth_monom, t), degree(kth_monom, u)
        kth_term_exps = filter(monom -> total_degree(monom) == t_deg, collect(monomials(Gk)))
        kth_term_exps_evaluated = map(monom -> evaluate(monom, point), kth_term_exps)
        kth_term_evaluated = kth_term_evaluated[1:length(kth_term_exps_evaluated)]
        @warn "" t_deg kth_term_exps kth_term_evaluated kth_term_exps_evaluated
        kth_term_coeffs = solve_transposed_vandermonde(
            ring_univariate,
            kth_term_exps_evaluated,
            kth_term_evaluated
        )
        interpolated_sum_of_terms = sum(kth_term_coeffs .* kth_term_exps)
        @warn "" interpolated_sum_of_terms
        interpolated_sum_of_terms = evaluate(interpolated_sum_of_terms, original_vars[1:k])
        interpolated_sum_of_terms = interpolated_sum_of_terms * last_var^(u_deg)
        Gn += interpolated_sum_of_terms
    end
    @info """
    The number of terms in the factors:
    Gn, Pn, Qn: $(length(Gn)), $(length(P)), $(length(Q))
    Points used: $T"""
    # push!(trace_info, (nvars=n, len_A=length(An), len_B=length(Bn), len_F=length(F)))
    @info "Interpolated" Gn
    Gn
end

####################
####################

# using Nemo

# R, (x1, x2, x3) = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:3]...], ordering=:degrevlex)

# P = ((x1 + 1) * (x1 + x2 * x3 + 1));
# Q = ((x1 + 2) * (x1 + x2 * x3 + 1));
# # P = ((x1 + 1) * (x1 * x2 * x3 + 1));
# # Q = ((x1 + 2) * (x1 * x2 * x3 + 1));

# G = iterative_gcd(P, Q);
# G_true = Nemo.gcd(P, Q);

# G
# G_true

# G == G_true
