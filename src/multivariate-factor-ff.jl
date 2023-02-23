
# Select the variable to be used as the main variable
# in the multivariate factorization.
# Currently, this defaults to the first variable.
#
# Returns a tuple 
# (index of main variable, main variable, all other variables)
function select_main_variable(F)
    R = parent(F)
    xs = gens(R)
    K = base_ring(R)
    for i in 1:length(xs)
        degree(F, xs[i]) < 2 && continue
        p = append!(zeros(R, i - 1), [xs[i]], zeros(R, length(xs) - i))
        f = evaluate(F, p)
        lc = leading_coefficient_in(F, xs[i])
        Runiv, _ = K["u"]
        funiv = to_univariate(Runiv, f)
        if !iszero(funiv) && !iszero(evaluate(lc, p))
            if degree(funiv) > 0
                if issquarefree(funiv)
                    return i, xs[i], append!(xs[1:i-1], xs[i+1:end])
                end
            end
        end
    end
    @warn "No variable provides a squarefree univariate factorization. Returning.."
    0, xs[1], xs[2:end]
end

# Given a factorization 
#   F(t,u) = lc(u)*P(t,u)*Q(t,u) 
# normalizes and returns the factors according to some consistent convention. 
#
# This currently assumes that F is monic as a polynomial in K[u][t]
# and normalizes P to be monic in K[u][t].
function normalize_factorization(Fi::Vector{T}) where {T}
    ai = map(trailing_coefficient, Fi)
    Fi = map(divexact, Fi, ai)
    Fi[1] = Fi[1]*prod(ai)
    Fi[1] = divexact(Fi[1], trailing_coefficient(Fi[1]))
    Fi
end

function costruct_evalpoint(t_idx, t, u, cs)
    n = length(cs) + 1
    @assert base_ring(t) == base_ring(u)
    @assert 1 <= t_idx <= n
    point = Vector{typeof(u)}(undef, n)
    point[t_idx] = t
    for i in 1:t_idx-1
        point[i] = u*cs[i]
    end
    for i in t_idx+1:n
        point[i] = u*cs[i-1]
    end
    point
end

# Returns some factors the given multivariate polynomial F
# over some prescibed finite field.
#
# T is the bound on the number of terms in the factors.
function multivariate_split_ff(F, T::Integer)
    @assert T > 0
    # P is a polynomial in K[xs]
    R = parent(F)
    K, xs = base_ring(R), gens(R)
    n = length(xs)
    @assert characteristic(K) > 0
    # Select the main variable
    t_idx, t, othervars = select_main_variable(F)
    # no good variable was found
    if t_idx == 0
        return [F]
    end
    @info "" t_idx t othervars
    # R_coeff is K[other].
    # R_u is K[t, u]
    R_coeff, _ = PolynomialRing(K, map(string, othervars))
    R_u, (t_u, u_u) = PolynomialRing(K, ["t", "u"])
    # Initialize the interpolation routine in the ring R_coeff
    # (it will interpolate the coefficients in other)    
    partial_degrees = [Nemo.degree(F, i) for i in 1:n if i != t_idx]
    # @info "" partial_degrees
    interpolator = FasterMultivariateBenOrTiwari(R_coeff, T, partial_degrees)
    # cs = [c2...cn] is an element of K^(n-1), which will be the starting
    # point of the geometric sequence as requested by the interpolator
    cs = startingpoint(interpolator)
    # Reduce P to a bivariate polynomial:
    # evaluate P at (t, c1^0*u, c2^0*u, ..., cn^0*u)
    # (the first point in the sequence will be c^0).
    point = costruct_evalpoint(t_idx, t_u, u_u, cs .^ 0)
    # @info "" cs point
    F_c_f = evaluate(F, point)
    # @assert isone(denominator(F_c_f))
    F_c = F_c_f
    # @info "Evaluated" F_c typeof(F_c)
    # Factor bivariate F_c.
    # (!) For now, we assume that there are only two factors.
    #
    # Reveals the bivariate factors P_c and Q_c, together with
    # univariate factors f1,...,fk and disjoint sets S_P and S_Q, so that
    #   F_c = P_c*Q_c
    #   F_c(t, u) = f1(t)*...*fk(t)  mod(u) 
    #   P_c(t, u) = prod_i f_i(t)  mod(u), i in S_P
    #   Q_c(t, u) = prod_j f_j(t)  mod(u), j in S_Q
    Fi, fi, Si = revealing_bivariate_factorization_ff(F_c)
    # @info "Revealing factorization" Fi fi Si
    Fi = normalize_factorization(Fi)
    # @info "Normalized" Fi
    # The storages to collect the coefficients of factors 
    # for different evaluation points
    Fcoeffs = Vector{Vector{Vector{elem_type(K)}}}(undef, length(Fi))
    # ..initialize storages and fill the values for the first evaluation point.
    for i in 1:length(Fi)
        Fcoeffs[i] = Vector{Vector{elem_type(K)}}(undef, length(Fi[i]))
        for j in 1:length(Fcoeffs[i])
            Fcoeffs[i][j] = Vector{elem_type(K)}(undef, 2T)
            Fcoeffs[i][j][1] = Nemo.coeff(Fi[i], j)
        end
    end
    # Evaluate another 2T - 1 times
    cpoints = map(i -> cs .^ i, 0:2T-1)
    for j in 2:2T
        point = costruct_evalpoint(t_idx, t_u, u_u, cpoints[j])
        F_c_f = evaluate(F, point)
        F_c = F_c_f
        # Here, all of these factorizations should benefit from the fact
        # that we already know correct in univariate factorization
        # and the corresponding recombination pattern.
        # We just need to hensel lift to F_c(t, u) mod <u^D>
        Fi = bivariate_factorization_ff_with_known_univariate(F_c, fi, Si)
        # @debug "Factor $j" Fi
        Fi = normalize_factorization(Fi)
        # @debug "Normalized $j" Fi
        @assert length(Fi) == length(Fcoeffs)
        # Collect the coefficients
        for i in 1:length(Fi)
            @assert length(Fi[i]) == length(Fcoeffs[i])
            for k in 1:length(Fcoeffs[i])
                Fcoeffs[i][k][j] = Nemo.coeff(Fi[i], k)
            end
        end
    end
    # @debug "Data collected" Fcoeffs
    # Interpolate !
    coeffs_interpolated = Vector{Vector{elem_type(R_coeff)}}(undef, length(Fcoeffs))
    for i in 1:length(Fcoeffs)
        coeffs_interpolated[i] = Vector{elem_type(R_coeff)}(undef, length(Fcoeffs[i]))
        for j in 1:length(Fcoeffs[i])
            coeffs_interpolated[i][j] = interpolate!(interpolator, cpoints, Fcoeffs[i][j])
        end
    end
    # @info "Interpolated coeffs" coeffs_interpolated
    # Collect the interpolated coefficients into a polynomial.
    Fmonoms = map(collect ∘ monomials, Fi)
    Finterpolated = map(
        (ms, cfs) -> 
            sum(evaluate(tt, [t, one(t)])*evaluate(cc, othervars) 
            for (tt, cc) in zip(ms, cfs)), 
        Fmonoms,
        coeffs_interpolated
    )
    Finterpolated = map(f -> evaluate(f, xs), Finterpolated)
    Finterpolated = map(f -> divexact(f, leading_coefficient(f)), Finterpolated)
    Finterpolated
end

# Returns some factors the given multivariate polynomial F
# over some prescibed finite field.
#
# T is the bound on the number of terms in the factors.
function multivariate_factor(F, T::Integer)
    @assert T > 0
    Fi = multivariate_split_ff(F, T)
    # @warn "" Fi
    length(Fi) == 1 && return Fi
    F2 = divexact(F, prod(Fi))
    isone(F2) && return Fi
    subfactors = multivariate_factor(F2, T) 
    union(Fi, subfactors)
end

#=

using Nemo

K = GF(2^31-1)
R, (x1,x2,x3,x4,x5,x6,x7) = PolynomialRing(K, ["x$i" for i in 1:7])

f = prod(x1 + x2^rand(0:2) + x3^rand(0:2) + rand(K) for _ in 1:3)*prod(x2^2 + x3^rand(0:2) + x4^rand(0:2) + rand(K) for _ in 1:3)*prod(x4 + x5^rand(0:2) + x6^rand(0:2) + x7^rand(0:2) + rand(K) for _ in 1:3)

=#

