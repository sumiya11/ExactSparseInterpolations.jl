
# Given F in K[x1,...,xn], returns F(x1,...,xk,ck+1,...,cn)
function project(F, k, c)
    ring = parent(F)
    @assert k < nvars(ring)
    _, new_vars = PolynomialRing(ring, gens(ring)[1:k])
    eval_point = vcat(new_vars, c[k+1:nvars(ring)])
    evaluate(F, eval_point)
end

function iterative_factor(F)
    ring = parent(F)
    field = base_ring(ring)
    c = distinct_points(field, n)
    F2 = project(F, 2, c)
    A2, B2 = Nemo.factor(F2)
    An, Bn = _iterative_factor(F, A2, B2, c)
    An, Bn
end

function _iterative_factor(F, A2, B2, c)
    ring = parent(F)
    field = base_ring(ring)
    n = nvars(ring)
    if n == 2
        return A2, B2
    end
    # Recursive call to obtain the factorization in n - 1 variables
    Fk = project(F, n - 1, c)
    Ak, Bk = _iterative_factor(Fk, A2, B2, c)
    # Evaluate F and Ak at a geometric progression
    T = length(F)
    interpolator = PrimesBenOrTiwari(ring, T)
    point = startingpoint(interpolator)
    F_coeffs_eval = simultaneous_multivariate_evaluate(ring, F_coeffs, point, 2T)
    P_coeffs_eval = simultaneous_multivariate_evaluate(ring, P_coeffs, point, 2T)
    # Hensel lift Pk and Qk to Pn and Qn
    for i in 1:2T
        Q_coeffs_eval = map(F_coeffs_eval)
        Pni = lift(P[i])
    end
    fast_vandermonde_solve(monoms, P_coeffs_eval)
end
