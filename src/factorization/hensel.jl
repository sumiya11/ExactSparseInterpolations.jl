# Hense lifting algorithms.

function leading_coefficient_in(f, x)
    @assert parent(f) == parent(x)
    lead_deg_in_x = Nemo.degree(f, x)
    idx_in_x = findfirst(==(x), gens(parent(f)))
    @assert !isnothing(idx_in_x)
    coeff(f, [idx_in_x], [lead_deg_in_x])
end

function content_in(f, x)
    Nemo.degree(f, x) == 0 && return one(f)
    @assert parent(f) == parent(x)
    idx_in_x = findfirst(==(x), gens(parent(f)))
    @assert !isnothing(idx_in_x)
    cont = coeff(f, [idx_in_x], [0])
    for i in 1:Nemo.degree(f, x)
        cont = gcd(cont, coeff(f, [idx_in_x], [i]))
    end
    cont
end

function primpart_in(f, x)
    cont = content_in(f, x)
    isone(cont) && return f
    divexact(f, cont * leading_coefficient(f))
end

function my_primpart(F)
    xs = vars(F)
    F_prim = F
    for x in xs
        F_prim = primpart_in(F_prim, x)
    end
    F_prim
end

# Given polynomials f, g, h ∈ K[y][x], and m ∈ K[y] such that
#   f - g h = 0 (mod m)
# and s g + t h = 1 (mod m),
# returns a quad (g*, h*, s*, t*) such that
#   f* - g* h* = 0 (mod m^2)
# and
#  s* g* + t* h* = 1 (mod m^2)
#
# # O( M(deg(f)) M(deg(m)) )
function hensel_step(f, g, h, s, t, m)
    m2 = m^2
    # well, this just calls bivariate polynomial division, which is sad
    e = mod(f - g * h, m2)
    q, r = divrem(s * e, h)
    gstar = mod(g + t * e + q * g, m2)
    hstar = mod(h + r, m2)
    b = mod(s * gstar + t * hstar - one(s), m2)
    c, d = divrem(s * b, hstar)
    sstar = mod(s - d, m2)
    tstar = mod(t - t * b - c * gstar, m2)
    gstar, hstar, sstar, tstar
end

# Given polynomials F, f1,...,fr ∈ K[y][x], p ∈ K[y], and integer l such that
#   F = lc_x(F)f1...fr    (mod p),
# returns polynomials (f1*, ..., fr*) such that
#   F* = lc_x(F)f1*...fr* (mod p^l)
#
# (!) Returned polynomials are monic.
# (!) Assumes p is linear, that is
#   p = y + a
# (!) Assumes variable `x` is the first variable in the polynomial ring object.
#
# O( log r (M(n)M(l) + M(n)log n) )
# where n = deg_x(F)
function hensel_multifactor_lifting(F::T, fi::Vector{T}, l::Integer, p::T) where {T}
    lcf = leading_coefficient_in(F, first(gens(parent(F))))
    @assert mod(lcf * prod(fi), p) == mod(F, p)
    r = length(fi)
    if isone(r)
        # this can be turned into inverting modulo p,
        # plus newton iteration
        lcfinv = invmod(lcf, p^l)
        return [mod(F * lcfinv, p^l)]
    end
    k = div(r, 2)
    d = ceil(Int, log2(l))
    gi = mod(lcf * prod(fi[i] for i in 1:k), p)
    hi = mod(prod(fi[i] for i in (k + 1):r), p)
    I_am_one, si, ti = Nemo.gcdx(gi, hi)
    @assert isone(I_am_one)
    si, ti = mod(si, p), mod(ti, p)
    for j in 1:d
        m = p^(2^(j - 1))
        gi, hi, si, ti = hensel_step(F, gi, hi, si, ti, m)
    end
    f1tok = hensel_multifactor_lifting(gi, fi[1:k], l, p)
    fktor = hensel_multifactor_lifting(hi, fi[(k + 1):r], l, p)
    append!(f1tok, fktor)
end

function hensel_normalization_joris(F::T, Fi::Vector{T}, l, p) where {T}
    R = parent(F)
    K = base_ring(R)
    x, y = gens(R)
    point = [R(rand(K)), y]
    R_uni, t = K["t"]
    result = similar(Fi)
    @debug """
    Lifted polys: $Fi
    Lead in x (x is main var.): $(leading_coefficient_in(F, x))
    Lifted up to $(p)^$(l)
    Evaluation point: $point"""
    for i in 1:length(Fi)
        fi = Fi[i]
        fi = evaluate(fi, point)
        fi_uni = to_univariate(R_uni, fi)
        r_, t_, s_ = Padé(fi_uni, t^(l), l - 2)
        B = evaluate(t_, y)
        B = divexact(B, leading_coefficient(B))
        result[i] = mod(Fi[i] * B, p^(l))
    end
    result
end

function hensel_lift(F::T, Fi::Vector{T}, l::Integer, p::T) where {T}
    R = parent(F)
    x, y = gens(R)
    @assert all(f -> isone(leading_coefficient_in(f, x)), Fi)
    Fi = hensel_multifactor_lifting(F, Fi, l, p)
    Fi = hensel_normalization_joris(F, Fi, l, p)
    Fi
end
