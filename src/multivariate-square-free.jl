# Multivariate square-free factorization over finite fields

# Returns the square-free factorization of P
function multivariate_square_free(F)
    ring = parent(F)
    field = base_ring(ring)
    # Random dilation, just in case
    dilation = distinct_points(field, nvars(ring))
    F = evaluate(F, dilation .* gens(ring))
    # Find a map from f(x_1,...,x_n) to f(x_1 t^i_1, ..., x_n t^i_n), such
    # that the images of P and Q are w-regular
    w = find_tagging_map(F)
    F_hat_coeffs = apply_tagging_map(F, w)
    @info "Tagging map: $w, weighted degree is $(length(F_hat_coeffs) - 1)"
    @debug "" F_hat_coeffs
    T = length(F)
    interpolator = PrimesBenOrTiwari(ring, T)
    point = startingpoint(interpolator)
    @info "Evaluation point is $point, the number of terms is $T"
    F_hat_coeffs_eval = simultaneous_multivariate_evaluate(ring, F_hat_coeffs, point, 2T)
    ring_univariate, t = PolynomialRing(field, "t")
    square_free_univariate_coeffs = Vector{Vector{Vector{elem_type(field)}}}(undef, 2T)
    lead = wlead(F, w)
    @info "The leading part is $lead"
    for i in 1:2T
        F_hat_eval = ring_univariate(F_hat_coeffs_eval[i])
        @info "" F_hat_eval
        square_free_hat = collect(Nemo.factor_squarefree(F_hat_eval))
        square_free_hat = map(collect, square_free_hat)
        @info "after factor" square_free_hat
        max_square_degree = maximum(last, square_free_hat)
        for j in 1:length(square_free_hat)
            square_free_hat[j][1] = shift_right(square_free_hat[j][1], valuation(square_free_hat[j][1], t))
            square_free_hat[j][1] = divexact(square_free_hat[j][1], leading_coefficient(square_free_hat[j][1]))
            @assert isone(leading_coefficient(square_free_hat[1][1]))
            @assert iszero(valuation(square_free_hat[j][1], t))
        end
        square_free_univariate_coeffs[i] = Vector{Vector{elem_type(field)}}(undef, max_square_degree)
        for j in 1:max_square_degree
            square_free_univariate_coeffs[i][j] = [one(field)]
        end
        @info "" square_free_hat
        for j in 1:length(square_free_hat)
            lead_monom_eval = evaluate(lead, point .^ (i - 1))
            factor_coeffs = collect(coefficients(square_free_hat[j][1]))
            # Correction: multiply by the evalution of a large monomial to keep
            # things polynomial.
            factor_coeffs .*= lead_monom_eval
            square_free_univariate_coeffs[i][square_free_hat[j][2]] = factor_coeffs
        end
    end
    # Interpolate each coefficient of the factorization
    @info "" square_free_univariate_coeffs
    nsquares = length(first(square_free_univariate_coeffs))
    multivariate_coeffs = Vector{Vector{elem_type(ring)}}(undef, nsquares)
    points = map(i -> point .^ i, 0:2T-1)
    factors_hat = Vector{elem_type(ring)}(undef, nsquares)
    for i in 1:length(multivariate_coeffs)
        ncoeffs = length(first(square_free_univariate_coeffs)[i])
        multivariate_coeffs[i] = Vector{elem_type(ring)}(undef, ncoeffs)
        for j in 1:length(multivariate_coeffs[i])
            factor_coeffs_ij = [cfs[i][j] for cfs in square_free_univariate_coeffs]
            multivariate_coeffs[i][j] = interpolate!(interpolator, points, factor_coeffs_ij)
        end
        @debug "" multivariate_coeffs[i]
        factors_hat[i] = sum(multivariate_coeffs[i])
    end
    @info "" factors_hat
    # Account for the valuation in x_1, ..., x_n
    xs = gens(ring)
    for i in 1:length(factors_hat)
        factor_hat = factors_hat[i]
        for (j, x) in enumerate(xs)
            s = i == valuation(F, xs[j])
            correction_degree = s - valuation(factors_hat[i], xs[j])
            if correction_degree < 0
                factor_hat = divexact(factor_hat, x^(-correction_degree))
            else
                factor_hat = factor_hat * x^(correction_degree)
            end
        end
        factors_hat[i] = factor_hat
    end
    true_count = findlast(factor -> !isone(factor), factors_hat)
    factors_hat = factors_hat[1:true_count]
    # Reverse dilation
    square_free_factors = map(factor -> evaluate(factor, inv.(dilation) .* gens(ring)), factors_hat)
    square_free_factors = map(factor -> divexact(factor, leading_coefficient(factor)), square_free_factors)
    square_free_factors
end

####################
####################

# For example, you can run the following commands in Julia.

# using Nemo

R, (x, y, z) = PolynomialRing(GF(2^62 + 135), ["x", "y", "z"])

F = x * y^3 * (y + 1)^3
Nemo.factor(F)

multivariate_square_free(F)

R, x = PolynomialRing(GF(2^62 + 135), [["x$i" for i in 1:6]...])

F = 6 * x[1] * x[5]^4 * x[6]^5 * (x[1] + 5x[6])^5

@time Nemo.factor(F)

@time sq = multivariate_square_free(F)
