
# Recombines the univariate factors Fi of bivariate F
# modulo the given modulo.
# (!) Assumes there are only two factors.
function recombine_factors(F::P, fi::Vector{P}, modulo::P) where {P}
    t, u = gens(parent(F))
    b = leading_coefficient_in(F, t)
    T = collect(1:length(fi))
    Fstar = F
    Fi = Vector{typeof(F)}()
    Si = Vector{Vector{Int}}()
    s = 1
    while 2s <= length(T)
        for S in Combinatorics.combinations(T, s)
            gstar = mod(b*prod(fi[i] for i in S), modulo)
            hstar = mod(b*prod(fi[i] for i in setdiff(T, S)), modulo)
            if total_degree(gstar * hstar) == total_degree(b*Fstar)
                T = setdiff(T, S)
                push!(Fi, primpart_in(gstar, t))
                push!(Si, S)
                Fstar = primpart_in(hstar, t)
                b = leading_coefficient_in(Fstar, t)
                s -= 1
                break
            end
        end
        s += 1
    end
    push!(Fi, Fstar)
    push!(Si, T)
    Fi, Si
end

# Computes a factorization of f(x, y) over a finite field,
# 
# (!) Assumes f(x, y) is monic & squarefree as a polynomial in K[y][x].
# (!) Assumes f(x, 0) is squarefree and lead_x(f)(y=0) is not zero.
# (!) assumes there are exactly two factors, P(x,y) and Q(x,y),
#     and that both factors depend on x.
function revealing_bivariate_factorization_ff(F::T, at_u=zero(base_ring(F))) where {T}
    R = parent(F)
    K, (t, u) = base_ring(R), gens(R)
    n, d = Nemo.degree(F, t), Nemo.degree(F, u)
    b = leading_coefficient_in(F, t)
    # modulo for reduction to univariate case
    m = u - at_u
    f = mod(F, m)
    @assert n >= 1
    @assert isone(gcd(f, derivative(f, t)))
    @assert !iszero(mod(b, m))
    Runiv, _ = K["t"] 
    # l is the bound on the length of hensel iteration
    l = d + 1 + Nemo.degree(b, u)
    # factor F mod m
    # into factors f1..fr
    funiv = to_univariate(Runiv, f)
    factorization = Nemo.factor(funiv)
    fiuniv = collect(keys(factorization.fac))
    # lift the factors f1..fr via Hensel lifting
    fi = map(f -> evaluate(f, t), fiuniv)
    filifted = hensel_multifactor_lifting(F, fi, l, m)
    # recombine the factors 
    # into true factors F1..Fk
    Fi, Si = recombine_factors(F, filifted, m^l)
    fi = map(S -> prod(fi[i] for i in S), Si)
    Fi, fi
end

function reveal_univariate_factors_ff(F::T, at_u=zero(base_ring(F))) where {T}
    R = parent(F)
    K, (t, u) = base_ring(R), gens(R)
    n, d = Nemo.degree(F, t), Nemo.degree(F, u)
    b = leading_coefficient_in(F, t)
    # modulo for reduction to univariate case
    m = u - at_u
    f = mod(F, m)
    @assert n >= 1
    @assert isone(gcd(f, derivative(f, t)))
    @assert !iszero(mod(b, m))
    Runiv, _ = K["t"] 
    # l is the bound on the length of hensel iteration
    l = d + 1 + Nemo.degree(b, u)
    # factor F mod m
    # into factors f1..fr
    funiv = to_univariate(Runiv, f)
    factorization = Nemo.factor(funiv)
    fiuniv = collect(keys(factorization.fac))
    fiuniv
end

function bivariate_factorization_ff_with_known_univariate(F, fi, at_u=zero(base_ring(F)))
    t, u = gens(parent(F))
    b = leading_coefficient_in(F, t)
    m = u - at_u
    n, d = Nemo.degree(F, t), Nemo.degree(F, u)
    l = d + 1 + Nemo.degree(b, u)
    # apply Hensel lifting
    Fi = hensel_multifactor_lifting(F, fi, l, m)
    Fi = map(f -> primpart_in(mod(b*f, m^l), t), Fi)
    Fi
end
