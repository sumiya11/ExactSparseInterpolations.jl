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

function multivariate_root_extraction(F)
    orig_ring = parent(F)
    ring = parent(F)
    field = base_ring(ring)
    # Random dilation, just in case
    dilation = distinct_points(field, nvars(ring))
    F = evaluate(F, dilation .* gens(ring))
    # Find a map from f(x_1,...,x_n) to f(x_1 t^i_1, ..., x_n t^i_n), such
    # that the images of P and Q are w-regular
    w = find_regularizing_weight(F)
    F_hat_coeffs = apply_regularizing_weight(F, w)
    @info """
    Tagging map: $w
    Initial total degree: $(max(total_degree(F)))
    Weighted degree: $(length(F_hat_coeffs))"""
    @assert length(F_hat_coeffs[end]) == 1
    @assert sum(F_hat_coeffs) == F
    # Evaluate F_hat at a geometric progression. 
    # Use the min. number of terms in the input as a crude estimation.
    T = length(F)
    interpolator = PrimesBenOrTiwari(ring, T)
    point = startingpoint(interpolator)
    @info """
    Evaluation point: $point, dilation: $dilation
    Interpolating for $T terms"""
    F_hat_coeffs_eval = simultaneous_multivariate_evaluate(ring, F_hat_coeffs, point, 2T)
    # Compute the root of the specialized F_hat
    ring_univariate, t = PolynomialRing(field, "t")
    roots_univariate_coeffs = Vector{Vector{elem_type(field)}}(undef, 2T)
    lead = wlead(F, w)
    for i in 1:(2T)
        F_hat_eval = ring_univariate(F_hat_coeffs_eval[i])
        root_hat = Nemo.sqrt(F_hat_eval)
        root_hat = shift_right(root_hat, valuation(root_hat, t))
        root_hat = divexact(root_hat, leading_coefficient(root_hat))
        @assert isone(leading_coefficient(root_hat))
        @assert iszero(valuation(root_hat, t))
        roots_univariate_coeffs[i] = collect(coefficients(root_hat))
        # Correction: multiply by the evalution of a large monomial to keep
        # things polynomial.
        lead_monom_eval = evaluate(lead, point .^ (i - 1))
        roots_univariate_coeffs[i] .*= lead_monom_eval
    end
    # Interpolate each coefficient of the gcd
    R_coeffs = Vector{elem_type(ring)}(undef, length(first(roots_univariate_coeffs)))
    points = map(i -> point .^ i, 0:(2T - 1))
    for i in 1:length(R_coeffs)
        root_coeffs_i = [cfs[i] for cfs in roots_univariate_coeffs]
        R_coeffs[i] = interpolate!(interpolator, points, root_coeffs_i)
    end
    R_hat = sum(R_coeffs)
    # Account for the valuation in x_1, ..., x_n
    xs = gens(ring)
    for x in xs
        correction_degree = div(valuation(F, x), 2) - valuation(R_hat, x)
        if correction_degree < 0
            R_hat = divexact(R_hat, x^(-correction_degree))
        else
            R_hat = R_hat * x^(correction_degree)
        end
    end
    R = R_hat
    # Reverse dilation
    R = evaluate(R, inv.(dilation) .* gens(ring))
    R = divexact(R, leading_coefficient(R))
    R
end

####################
####################

# # using Nemo, ExactSparseInterpolations
# using Nemo

# R, (x1, x2, x3, x4, x5, x6) = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:6]...])

# F = (x1 * x2 - 1)^2;
# length(F), degrees(F)

# @time R1 = Nemo.sqrt(F);
# @time R2 = ExactSparseInterpolations.multivariate_root_extraction(F);

# R1
# R2

# length(R1), degrees(R1)
# R1 == R2
