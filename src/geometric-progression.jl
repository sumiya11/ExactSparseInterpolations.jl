
function fastmultivariateinterpolate(R, xs, ys)
    @assert length(xs) == length(ys)
    K = base_ring(R)
    d = length(xs)
    # ω = distinctrootsofunity(K, d)
    _, L, _ = fastconstrainedEEA()
    mi = Nemo.roots(Lrev)
    # find the monomials of the interpolant:
    monoms = [
        discrete_log(m, bot.ω0)
        for m in mi
    ]
    # find the coefficients of the interpolant:
    # matrix A is a t×t Vandermonde matrix,
    # solve A*ai = vi for ai
    Vspace = Nemo.MatrixSpace(base_ring(R), t, t)
    A = Vspace(reshape([
        mi[i] ^ (j)
        for i in 1:t
            for j in 1:t
    ], (t, t)))
    vi = bot.vi[1:t]
    ai = Nemo.inv(A) * vi
    true, R(ai, monoms)
end
