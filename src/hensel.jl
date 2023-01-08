
function hensel_step(f, g, h, s, t, m)
    m2 = m^2
    e = mod(f - g*h, m2)
    q, r = divrem(s*e, h)   # should be reduced mod m^2
    gstar = mod(g + t*e + q*g, m2)
    hstar = mod(h + r, m2)
    b = mod(s*gstar + t*hstar - one(s), m2)
    c, d = divrem(s*b, hstar)  # should be reduced mod m^2
    sstar = mod(s - d, m2)
    tstar = mod(t - t*b - c*gstar, m2)
    gstar, hstar, sstar, tstar
end

function hensel_multifactor_lifting(f, fi, l, p)
    lcf = leading_coefficient(f)
    @assert mod(lcf*prod(fi), p) == mod(f, p)
    r = length(fi)
    isone(r) && return [mod(f, p^l)]
    k = div(r, 2)
    d = ceil(Int, log2(l))
    gi = mod(lcf*prod(fi[i] for i in 1:k), p)
    hi = mod(lcf*prod(fi[i] for i in k+1:r), p)
    _, si, ti = Nemo.gcdx(gi, hi)
    si = mod(si, p)
    ti = mod(ti, p)
    for j in 1:d
        m = p^(2^(j-1))
        gi, hi, si, ti = hensel_step(f, gi, hi, si, ti, m)
    end
    f1tok = hensel_multifactor_lifting(gi, fi[1:k], l, p)
    fktor = hensel_multifactor_lifting(hi, fi[k+1:r], l, p)
    append!(f1tok, fktor)
end

# f = 1 * (x + y) * (x^2 + 2*y) * (x + 3*y^2)

function _bivariate_factorization_ff(f)
    R = parent(f)
    K, (x, y) = base_ring(R), gens(R)
    n, d = degree(f, x), degree(f, y)
    @assert isone(gcd(f, derivative(f, x)))
    isone(n) && return [f]
    @assert order(K) >= 4n*d
    Runiv, _ = K["u"]
    b = leading_coefficient(f)  # coefficient in x
    # @assert isconstant(b)
    u = rand(K)
    while true
        u = rand(K)
        fbar = evaluate(f, [x, u])  # mod y - u
        iszero(b) && continue
        @info "" Nemo.factor(f)
        @info "" Nemo.factor(fbar) Nemo.factor(derivative(fbar, x))
        !isone(gcd(fbar, derivative(fbar, x))) && continue
        break
    end
    l = d + 1 + 0  # 0 == degree of b in y
    funiv = to_univariate(Runiv, mod(f, y - u))
    @info "" u f Nemo.factor(f) Nemo.factor(funiv)
    hi = map(kv -> kv[1]^kv[2], collect(factorization_ff(funiv)))
    hi = map(h -> evaluate(h, x), hi)
    gi = hensel_multifactor_lifting(f, hi, l, y - u)
    T = collect(1:length(gi))
    G = elem_type(R)[]
    fstar = f
    s = 1
    ml = (y - u)^l
    while 2s <= length(T)
        for S in Combinatorics.combinations(T, s)
            gstar = mod(b*prod(gi[i] for i in S), ml)
            hstar = mod(b*prod(gi[i] for i in setdiff(T, S)), ml)
            @info "" gstar hstar fstar
            if degree(gstar * hstar, x) == degree(b*fstar, x)
                T = setdiff(T, S)
                push!(G, gstar) # primpart(gstar))
                fstar = hstar # primpart(hstar)
                b = leading_coefficient(fstar)
                s -= 1
                break
            end
        end
        s += 1
    end
    push!(G, fstar)
end

function bivariate_factorization_ff(f)
    
end

