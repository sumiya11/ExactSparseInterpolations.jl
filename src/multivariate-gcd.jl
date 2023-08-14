# Multivariate gcd over finite fields

default_weight(n) = ones(Int, n)
wdegree(F) = wdegree(F, default_weight(nvars(parent(F))))
wdegree(F, w) = maximum(x -> dot(x, w), collect(exponent_vectors(F)))
wvaluation(F) = wvaluation(F, default_weight(nvars(parent(F))))
wvaluation(F, w) = minimum(x -> dot(x, w), collect(exponent_vectors(F)))

wlead(F) = wlead(F, default_weight(nvars(parent(F))))
wlead(F, w) = monomial(F, findmax(monom -> wdegree(monom, w), collect(monomials(F)))[2])

function homogenize(F)
    ring = parent(F)
    newstrings = vcat("x0", map(string, gens(ring)))
    newring, _ = PolynomialRing(base_ring(ring), newstrings)
    deg = total_degree(F)
    newmonoms = map(e -> vcat(deg - sum(e), e), collect(exponent_vectors(F)))
    newcoeffs = collect(coefficients(F))
    newring(newcoeffs, newmonoms)
end

# Returns a vector of regularizing weights
function find_tagging_map(P, Q)
    # An exponent vector with the largest norm in P and Q
    wP = find_tagging_map(P)
    wQ = find_tagging_map(Q)
    if wdegree(P, wP) > wdegree(P, wQ)
        w = wP
    else
        w = wQ
    end
    w
end

function find_tagging_map(P)
    normP, idxP = findmax(norm, collect(exponent_vectors(P)))
    w = exponent_vector(P, idxP)
    w
end

# Returns a vector of coefficients of the densification of F in w:
# [F_0, F_1, ..., F_d],
# where F_i is the sum of all terms of F with w-degree equal to i.
function apply_tagging_map(F, w)
    exps = collect(exponent_vectors(F))
    wdeg = wdegree(F, w)
    multivariate_coeffs = map(_ -> zero(F), 1:(wdeg+1))
    for i in 1:length(exps)
        wdeg_i = dot(exps[i], w)
        multivariate_coeffs[wdeg_i+1] += term(F, i)
    end
    multivariate_coeffs
end

# Returns the gcd of P and Q
function multivariate_gcd(P, Q; do_homogenize=false)
    orig_ring = parent(P)
    if do_homogenize
        d = max(total_degree(P), total_degree(Q))
        P, Q = homogenize(P), homogenize(Q)
        @info """
        Degree before homogenization: $d
        Degree after homogenization: $(max(total_degree(P), total_degree(Q)))"""
    end
    ring = parent(P)
    field = base_ring(ring)
    # Random dilation, just in case
    dilation = distinct_points(field, nvars(ring))
    P = evaluate(P, dilation .* gens(ring))
    Q = evaluate(Q, dilation .* gens(ring))
    # Find a map from f(x_1,...,x_n) to f(x_1 t^i_1, ..., x_n t^i_n), such
    # that the images of P and Q are w-regular
    w = find_tagging_map(P, Q)
    P_hat_coeffs = apply_tagging_map(P, w)
    Q_hat_coeffs = apply_tagging_map(Q, w)
    @info """
    Tagging map: $w
    Initial total degree: $(max(total_degree(P), total_degree(Q)))
    Weighted degree: $(max(length(P_hat_coeffs), length(Q_hat_coeffs)) - 1)"""
    @info "" P_hat_coeffs Q_hat_coeffs
    @assert length(P_hat_coeffs[end]) == 1 || length(Q_hat_coeffs[end]) == 1
    @assert sum(P_hat_coeffs) == P && sum(Q_hat_coeffs) == Q
    # Evaluate P_hat and Q_hat at a geometric progression. 
    # Use the min. number of terms in the input as a crude estimation.
    T = min(length(P), length(Q))
    interpolator = PrimesBenOrTiwari(ring, T)
    point = startingpoint(interpolator)
    @info """
    Evaluation point: $point, dilation: $dilation
    Interpolating for $T terms"""
    P_hat_coeffs_eval = simultaneous_multivariate_evaluate(ring, P_hat_coeffs, point, 2T)
    Q_hat_coeffs_eval = simultaneous_multivariate_evaluate(ring, Q_hat_coeffs, point, 2T)
    # Compute the gcd of the specialized P_hat and Q_hat
    ring_univariate, t = PolynomialRing(field, "t")
    gcds_univariate_coeffs = Vector{Vector{elem_type(field)}}(undef, 2T)
    lead = wlead(P, w) * wlead(Q, w)
    for i in 1:2T
        P_hat_eval = ring_univariate(P_hat_coeffs_eval[i])
        Q_hat_eval = ring_univariate(Q_hat_coeffs_eval[i])
        gcd_hat = gcd(P_hat_eval, Q_hat_eval)
        gcd_hat = shift_right(gcd_hat, valuation(gcd_hat, t))
        @assert isone(leading_coefficient(gcd_hat))
        @assert iszero(valuation(gcd_hat, t))
        gcds_univariate_coeffs[i] = collect(coefficients(gcd_hat))
        # Correction: multiply by the evalution of a large monomial to keep
        # things polynomial.
        # NOTE: this was not in the manuscript
        lead_monom_eval = evaluate(lead, point .^ (i - 1))
        gcds_univariate_coeffs[i] .*= lead_monom_eval
    end
    # Interpolate each coefficient of the gcd
    G_coeffs = Vector{elem_type(ring)}(undef, length(first(gcds_univariate_coeffs)))
    points = map(i -> point .^ i, 0:2T-1)
    for i in 1:length(G_coeffs)
        gcd_coeffs_i = [cfs[i] for cfs in gcds_univariate_coeffs]
        G_coeffs[i] = interpolate!(interpolator, points, gcd_coeffs_i)
    end
    G_hat = sum(G_coeffs)
    # Account for the valuation in x_1, ..., x_n
    xs = gens(ring)
    for x in xs
        correction_degree = min(valuation(P, x), valuation(Q, x)) - valuation(G_hat, x)
        if correction_degree < 0
            G_hat = divexact(G_hat, x^(-correction_degree))
        else
            G_hat = G_hat * x^(correction_degree)
        end
    end
    G = G_hat
    # Reverse dilation
    G = evaluate(G, inv.(dilation) .* gens(ring))
    if do_homogenize
        G = evaluate(G, vcat(one(orig_ring), gens(orig_ring)))
    end
    G = divexact(G, leading_coefficient(G))
    G
end

####################
####################

# For example, you can run the following commands in Julia.



# # First, install / update the packages:
# using Pkg
# Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
# Pkg.add("Nemo")

# Then, you should be able to do:
# using Nemo, ExactSparseInterpolations
# using Nemo, Combinatorics

# R, (x1, x2, x3, x4) = PolynomialRing(GF(2^62 + 135), ["x1", "x2", "x3", "x4"])

# random_monomial(xs, d) = prod(rand(xs) for i in 1:d)

# begin
#     for i in 1:100
#         xx = [x1, x2, x3]
#         D = rand(1:2)
#         sq1, sq2, sq3 = rand(1:1), rand(1:1), rand(1:1)
#         T1, T2, T3 = rand(1:10), rand(1:10), rand(1:10)
#         monoms1 = map(_ -> random_monomial(xx, D), 1:T1)
#         coeffs1 = map(_ -> rand(base_ring(R)), 1:length(monoms1))
#         monoms2 = map(_ -> random_monomial(xx, D), 1:T3)
#         coeffs2 = map(_ -> rand(base_ring(R)), 1:length(monoms2))
#         monoms3 = map(_ -> random_monomial(xx, D), 1:T3)
#         coeffs3 = map(_ -> rand(base_ring(R)), 1:length(monoms3))
#         poly1 = sum(coeffs1 .* monoms1)^sq1
#         poly2 = sum(coeffs2 .* monoms2)^sq2
#         poly3 = sum(coeffs3 .* monoms3)^sq3
#         # w = _find_tagging_map_hom(poly3 * poly1, poly3 * poly2)
#         w = _find_tagging_map_hom(poly3 * poly1)
#         P_hat = apply_tagging_map(poly3 * poly1, w)
#         # Q_hat = apply_tagging_map(poly3 * poly2, w)
#         @info "" poly3 * poly1 P_hat[end]
#         @info "" w
#         @assert length(P_hat[end]) == 1
#         # G1 = multivariate_gcd(poly3 * poly1, poly3 * poly2, do_homogenize=false)
#         # G2 = Nemo.gcd(poly3 * poly1, poly3 * poly2)
#         # @info "" G1 G2
#         # @assert G1 == G2
#     end
# end

# F = x1 * x2 + x1 * x3

# w = _find_tagging_map_hom(F)
# apply_tagging_map(F, w)

# G = (x1 + x2 + 5)^3

# homogenize(G)

# G == ExactSparseInterpolations.multivariate_gcd(F, G)

# G

# w = _find_tagging_map_hom(F, G)
# apply_tagging_map(F, w)

# R, x = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:10]...])

# P = (x[1] + x[2] + 2) * sum(x)^2 * (prod(x) - 1) * x[9]
# Q = (x[1] + x[2] + 3) * sum(x)^3 * (prod(x) + 1) * x[9] * x[10]

# G1 = Nemo.gcd(P, Q)
# G2 = ExactSparseInterpolations.multivariate_gcd(P, Q)

# @assert G1 == G2

