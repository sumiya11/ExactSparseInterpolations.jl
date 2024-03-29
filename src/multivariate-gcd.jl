# Multivariate gcd over finite fields

default_weight(n) = ones(Int, n)
wdegree(F) = wdegree(F, default_weight(nvars(parent(F))))
wdegree(F, w) = maximum(x -> dot(x, w), collect(exponent_vectors(F)), init=-1)
wvaluation(F) = wvaluation(F, default_weight(nvars(parent(F))))
wvaluation(F, w) = minimum(x -> dot(x, w), collect(exponent_vectors(F)), init=Inf)

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
function find_regularizing_weight(P, Q)
    # An exponent vector with the largest norm in P and Q
    wP = find_regularizing_weight(P)
    wQ = find_regularizing_weight(Q)
    if wdegree(P, wP) > wdegree(P, wQ)
        w = wP
    else
        w = wQ
    end
    w
end

function find_regularizing_weight(P)
    normP, idxP = findmax(norm, collect(exponent_vectors(P)))
    w = exponent_vector(P, idxP)
    w
end

# Returns a vector of coefficients of the densification of F in w:
# [F_0, F_1, ..., F_d],
# where F_i is the sum of all terms of F with w-degree equal to i.
function apply_regularizing_weight(F, w)
    ring = parent(F)
    exps = collect(exponent_vectors(F))
    wdeg = wdegree(F, w)
    densified_coeffs = map(_ -> Vector{elem_type(base_ring(ring))}(), 1:(wdeg + 1))
    densified_exps = map(_ -> Vector{Vector{Int}}(), 1:(wdeg + 1))
    for i in 1:length(exps)
        wdeg_i = dot(exps[i], w)
        push!(densified_coeffs[wdeg_i + 1], coeff(F, i))
        push!(densified_exps[wdeg_i + 1], exponent_vector(F, i))
    end
    multivariate_polys = map(ce -> ring(ce[1], ce[2]), zip(densified_coeffs, densified_exps))
    multivariate_polys
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
    w = find_regularizing_weight(P, Q)
    P_hat_coeffs = apply_regularizing_weight(P, w)
    Q_hat_coeffs = apply_regularizing_weight(Q, w)
    @info """
    Tagging map: $w
    Initial total degree: $(max(total_degree(P), total_degree(Q)))
    Weighted degree: $(max(length(P_hat_coeffs), length(Q_hat_coeffs)) - 1)"""
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
    for i in 1:(2T)
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
    points = map(i -> point .^ i, 0:(2T - 1))
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
#
# # Then, you should be able to do:
# using Nemo, ExactSparseInterpolations
#
# R, x = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:6]...])
#
# P = (x[1] + x[2] + 2) * sum(x)^2 * (prod(x) - 1) * x[5]
# Q = (x[1] + x[2] + 3) * sum(x)^3 * (prod(x) + 1) * x[5] * x[6]
#
# G1 = Nemo.gcd(P, Q)
# G2 = ExactSparseInterpolations.multivariate_gcd(P, Q)
#
# @assert G1 == G2
