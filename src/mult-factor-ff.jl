
# Factors the given multivariate polynomial P 
# over some prescibed finite field.
#
# (!) Assumes P factors into exactly 2 irreducible polynomials,
# and that each of these factors depends on the first indeterminate.
# (!) Assumes P is monic in the first variable. 
#
# D is the maximum of partial degrees in the polynomial factors.
# ( D can be inferred from the given P as
#   max(deg(P, x1), deg(P, x2), ...)),
#   but I want to be explicit here )
# T is the bound on the number of terms in the factors.
# (Can T be inferred from the given P ?? -- no)
function _multivariate_factorization_ff(P, T::Integer, D::Integer)
    # I use Nemo.factor here to factor bivariate polynomials,
    # so that we do not benefit from multiple successive calls 
    # to factorization yet.
    #
    # We assume that the input is monic in the first variable,
    # but the "the normalization problem" can still occur here.
    # Indeed, if
    #   P = Q*R
    # as long as lc(Q)*lc(R) = 1 in K, Q and R are unique only up to a unit
    # 
    # P is a polynomial in K[xs]
    R = parent(P)
    K, xs = base_ring(R), gens(R)
    n = length(xs)
    @assert characteristic(K) > 0
    @assert isone(leading_coefficient(P))
    # Our main variable t is just the first variable.
    t = first(xs)
    # R_coeff is K[x2..xn].
    # R_u is K[t, u]
    R_coeff, _ = PolynomialRing(K, map(string, xs[2:end]))
    R_u, (t_u, u) = PolynomialRing(K, ["t", "u"])
    # Initialize the interpolation routine in the ring R_coeff
    # (it will interpolate the coefficients in x2,..,xn)    
    partial_degrees = [D for _ in 2:n]
    interpolator = FasterMultivariateBenOrTiwari(R_coeff, T, partial_degrees)
    # cs = [c1...cn] is an element of K^(n-1), which will be the starting
    # point of the geometric sequence as requested by the interpolator
    cs = startingpoint(interpolator)
    # Reduce P to a bivariate polynomial:
    # evaluate P at (t, c1^0*u, c2^0*u, ..., cn^0*u)
    # (the first point in the sequence will be c^0).
    # The following syntax means:
    #   cs .^ 0 := [c1^0, c2^0, ..., cn^0]
    #   u .* (cs .^ 0) := [c1^0*u, c2^0*u, ..., cn^0*u]
    #   point := [t_u, c1^0*u, c2^0*u, ..., cn^0*u]
    point = [t_u, (u .* (cs .^ 0))...]
    P_c = evaluate(P, point)
    # Factor P_c into Q_c and R_c.
    # (!) For now, we assume that there are only two factors
    Q_c, R_c = Nemo.factor(P_c)
    Q_c, R_c = first(Q_c), first(R_c)
    # We do not control the order of returned factors from Nemo, 
    # so here is a quick hack that we will apply to make sure 
    # (at least for factors of different length)
    # that Q_c and R_c are in the same order for different evaluation points
    Q_c, R_c = sort([Q_c, R_c], by = length)
    # The storages to collect the coefficients of factors 
    # for different evaluation points
    Q_coeffs = Vector{Vector{elem_type(K)}}(undef, length(Q_c))
    R_coeffs = Vector{Vector{elem_type(K)}}(undef, length(R_c))
    # ..initialize storages and fill the values for the first evaluation point.
    # `@inbounds` tells to ommit all array access inbounds checking
    for (factor, storage) in ((Q_c, Q_coeffs), (R_c, R_coeffs))
        @inbounds for i in 1:length(storage)
            storage[i] = Vector{elem_type(K)}(undef, 2T)
            storage[i][1] = Nemo.coeff(factor, i)
        end
    end
    # Evaluate another 2T - 1 times
    points = map(i -> cs .^ i, 0:2T-1)
    for j in 2:2T
        point = [t_u, (u .* points[j])...]
        P_c = evaluate(P, point)
        # Here, all of these factorizations should benefit from the fact
        # that we already know correct pattern in bivariate factorization,
        # as well as correct pattern in the corresponding univariate factorizations.
        # So, instead of doing all four steps:
        #   - specialize P_c(t, 0)
        #   - factor P_c(t, 0)
        #   - hensel lift to P_c(t, y) mod <y^D>
        #   - combine the factors of P_c(t, y) mod <y^D>
        # we just need the first three.
        Q_c, R_c = Nemo.factor(P_c)
        Q_c, R_c = first(Q_c), first(R_c)
        Q_c, R_c = sort([Q_c, R_c], by = length)
        # Collect the coefficients
        for (factor, storage) in ((Q_c, Q_coeffs), (R_c, R_coeffs))
            for i in 1:length(storage)
                storage[i][j] = Nemo.coeff(factor, i)
            end
        end
    end
    # Interpolate !
    Q_interpolated = Vector{elem_type(R_coeff)}(undef, length(Q_c))
    R_interpolated = Vector{elem_type(R_coeff)}(undef, length(R_c))
    for i in 1:length(Q_coeffs)
        Q_interpolated[i] = interpolate!(interpolator, points, Q_coeffs[i])        
    end
    for i in 1:length(R_coeffs)
        R_interpolated[i] = interpolate!(interpolator, points, R_coeffs[i])        
    end
    # Collect the interpolated coefficients into a polynomial.
    # Q = sum t^i * C_i,
    # where t is in K[t, u],
    # and C_i is in K[x2..xn].
    # This evaluation is just a literal substitution 
    # for compatibility between different rings.
    Q_monoms = collect(monomials(Q_c))
    R_monoms = collect(monomials(R_c))
    Q = sum(evaluate(qt, [t, one(t)])*evaluate(qc, xs[2:end]) 
            for (qt, qc) in zip(Q_monoms, Q_interpolated))
    R = sum(evaluate(rt, [t, one(t)])*evaluate(rc, xs[2:end]) 
            for (rt, rc) in zip(R_monoms, R_interpolated))
    @assert Q*R == P
    Q, R
end

#=
Example:

using Nemo, ExactSparseInterpolations

R, (x,y,z,w) = GF(2^16+1)["x","y","z","w"]

P = (x^2 + 2x*y^2*z^2 + y*z*w + 17z^8)*(x^4 + 3z^5*w^4)

Q, R = ExactSparseInterpolations._multivariate_factorization_ff(P, 4, 8)

## prints
(x^4 + 3*z^5*w^4, x^2 + 2*x*y^2*z^2 + y*z*w + 17*z^8)

=#
