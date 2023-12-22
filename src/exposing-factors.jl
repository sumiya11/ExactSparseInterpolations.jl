using LinearAlgebra

"""
    get_exposing_monomial_map(F)

Returns a map of monomials of the form 

    xi => xi * x1^Ni,

such that in the image of `F` under this map all of the non-constant polynomial
factors depend on `x1`.
"""
function get_exposing_monomial_map(F, monomial_map_type)
    if monomial_map_type === :deterministic
        get_exposing_monomial_map_deterministic(F)
    elseif monomial_map_type === :regularizing
        get_exposing_monomial_map_regularizing(F)
    elseif monomial_map_type === :stochastic
        get_exposing_monomial_map_stochastic(F)
    end
end

function get_exposing_monomial_map_deterministic(F)
    ring = parent(F)
    n, xs = nvars(ring), gens(ring)
    d = total_degree(F)
    new_ring, xs = LaurentPolynomialRing(base_ring(ring), map(string, xs))
    x1 = xs[1]
    monomial_map = Dict()
    monomial_map[x1] = x1
    for i in 2:n
        xi = xs[i]
        deg = d * Primes.nextprime(2, i - 1)
        monomial_map[xi] = xi * x1^deg
    end
    # @info """
    # Polynomial in $n variables of total degree $d.
    # Monomial map:\n$monomial_map"""
    new_ring, monomial_map
end

function get_exposing_monomial_map_regularizing(F)
    ring = parent(F)
    w = argmax(m -> norm(m), collect(exponent_vectors(F)))
    n, xs = nvars(ring), gens(ring)
    d = total_degree(F)
    new_ring, xs = LaurentPolynomialRing(base_ring(ring), map(string, xs))
    x1 = xs[1]
    monomial_map = Dict()
    monomial_map[x1] = x1
    for i in 2:n
        xi = xs[i]
        deg = 2w[i]
        monomial_map[xi] = xi * x1^deg
    end
    @info """
    Polynomial in $n variables of total degree $d.
    Monomial map:\n$monomial_map"""
    new_ring, monomial_map
end

function get_exposing_monomial_map_stochastic(F)
    ring = parent(F)
    n, xs = nvars(ring), gens(ring)
    new_ring, xs = LaurentPolynomialRing(base_ring(ring), map(string, xs))
    d = total_degree(F)
    x1 = xs[1]
    monomial_map = Dict()
    monomial_map[x1] = x1
    for i in 2:n
        xi = xs[i]
        # deg = rand(vcat(-5:-1:-1, 1:1:5))
        deg = rand([-1, 1]) * rand(1:2)
        monomial_map[xi] = xi * x1^deg
    end
    # if n > 1
    #     monomial_map[xs[2]] = xs[2] * x1^2
    # end
    # if n > 2
    #     monomial_map[xs[3]] = xs[3] * x1^4
    # end
    # @info """
    # Polynomial in $n variables of total degree $d.
    # Monomial map:\n$monomial_map"""
    new_ring, monomial_map
end

"""
    apply_monomial_map(F, monomial_map)

Applies the monomial map `monomial_map` of the form
    
    xi => xi * x1^di,

to `F` and returns the resulting polynomial.
"""
function apply_monomial_map(F, new_ring, monomial_map)
    new_poly = zero(new_ring)
    for (i, term) in enumerate(terms(F))
        new_term = one(new_ring) * Nemo.coeff(F, i)
        for var in vars(term)
            exp = degree(term, var)
            exp == 0 && continue
            var = evaluate(var, gens(new_ring))
            new_var = monomial_map[var]
            new_term *= new_var^exp
        end
        new_poly += new_term
    end
    new_poly
end

function apply_inverse_of_monomial_map(F, orig_ring, monomial_map)
    ring = parent(F)
    xs = gens(ring)
    orig_xs = gens(orig_ring)
    invmap = Dict(xs[1] => xs[1])
    for (xvar, varimage) in monomial_map
        if xvar == xs[1]
            continue
        end
        varimage_orig = evaluate(varimage, orig_xs)
        d = -degree(varimage_orig, orig_xs[1])
        invmap[xvar] = xvar * xs[1]^d
    end
    new_poly = zero(ring)
    F_ = evaluate(F, orig_xs)
    for (i, term) in enumerate(terms(F_))
        new_term = one(ring) * Nemo.coeff(F_, i)
        for var in vars(term)
            exp = degree(term, var)
            exp == 0 && continue
            var = evaluate(var, xs)
            new_var = invmap[var]
            new_term *= new_var^exp
        end
        new_poly += new_term
    end
    new_poly
end

function factor_via_exposing(F, monomial_map_type=:deterministic)
    total_degree(F) < 1 && return true
    ring = parent(F)
    xs = gens(ring)
    # Remove the valuation
    for x in xs
        _, F = remove(F, x)
    end
    # Get and apply a monomial map
    F_tilde = F # evaluate(F, vcat(one(ring), xs[2:end]))
    new_ring, monomial_map = get_exposing_monomial_map(F_tilde, monomial_map_type)
    F_hat = apply_monomial_map(F_tilde, new_ring, monomial_map)
    # for x in xs
    #     _, F_hat = remove(F_hat, x)
    # end
    # Factor!
    success = true
    factors = Nemo.factor(F_hat)
    D = 2 * total_degree(F)^2
    a1 = sort(
        map(
            s -> remove(s, xs[1])[2],
            map(
                s -> evaluate(
                    apply_inverse_of_monomial_map(xs[1]^D * s[1], ring, monomial_map)^s[2],
                    xs
                ),
                collect(factors)
            )
        ),
        by=leading_monomial
    )
    a2 = sort(map(s -> s[1]^s[2], collect(Nemo.factor(F))), by=leading_monomial)
    println(a1 == a2)
    for f in factors
        fact, multiplicity = f
        # Success is all of the factors depend on x1
        if degree(evaluate(fact, xs), xs[1]) < 1
            success = false
            println("======================")
            println(Nemo.factor(F))
            println(factors)
            println("======================")
        end
    end
    uwu = prod(
        map(s -> remove(evaluate(s[1], xs)^s[2], xs[1])[end], collect(factors));
        init=one(ring)
    )
    @info """
    Degrees
    Before map: $(degrees(F)), total: $(total_degree(F))
    After map: $(degrees(uwu)), total: $(total_degree(uwu))"""
    success
end

#############################################
#############################################

# using Nemo

# n = 16
# d = 2
# R, xi = Nemo.PolynomialRing(QQ, [["x$i" for i in 1:n]...])

# factor_via_exposing((xi[2] * xi[3] + 2), :regularizing)

# begin
#     counter = 0
#     boot = 1000

#     for i in 1:boot
#         P1 = rand(R, 3:5, 0:d, 0:10)
#         P2 = rand(R, 3:4, 0:(d - 1), 0:3)
#         P3 = rand(R, 3:4, 0:(d - 1), 0:3)
#         P4 = rand(R, 3:4, 0:(d - 1), 0:3)
#         P5 = rand(R, 3:4, 0:(d - 1), 0:3)
#         P6 = xi[3]^rand(1:2) * xi[4]^rand(0:2) + rand(1:5) * rand(xi) + rand(1:5)
#         u1 = xi[end - 2]^rand(0:3) - rand(1:10)
#         u2 = xi[end - 1]^rand(0:3) - rand(1:10)
#         F = P1 * P2 * P3 * P4 * P5 * P6 * u1 * u2

#         success = factor_via_exposing(F, :regularizing)
#         if !success
#             @warn "Beda!"
#             counter += 1
#         end
#     end

#     if counter > 0
#         @error "Some errors detected: $counter out of $boot"
#     else
#         @info "All is okay!"
#     end
# end
