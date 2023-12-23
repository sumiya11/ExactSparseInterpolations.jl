
"""
    irreducible_factorization(F)

Returns a tuple (`c`, `factors`), where `c` is a unit and `factors` are powers
of irreducible polynomials, and the following condition holds:

    F == c * prod(factors)

"""
function irreducible_factorization(F)
    # 0. Remove monomial multiples
    ring, ground_field = parent(F), base_ring(F)
    monomial_multiples = Vector{typeof(F)}()
    for x in gens(ring)
        v, F = remove(F, x)
        if !iszero(v)
            push!(monomial_multiples, x^v)
        end
    end
    isconstant(F) && return coeff(F, 1), monomial_multiples

    # 0. Apply a random dilation (or dilatation?..)
    dilation = map(_ -> rand(ground_field), 1:nvars(ring))
    F = evaluate(F, gens(ring) .* dilation)
    @info "Using random dilation $(gens(ring)) -> $(gens(ring) .* dilation)"

    # 1. Apply the transformation, together with a regularizing map x_i -> x_i ξ^w_n
    α = map(ground_field, Primes.nextprimes(2, nvars(ring)))
    β = map(_ -> rand(ground_field), 1:nvars(ring))
    γ = map(_ -> rand(ground_field), 1:nvars(ring))
    F̂_regular = homotopy_substitution_with_regularizing_transform(F, α, β, γ)
    @debug "After the substitution (factored for convenience):" Nemo.factor(F̂_regular)

    # 2. Compute a 4-variate factorization of the first specialization
    new_ring = parent(F̂_regular)
    (ξ, xs..., t, u, v, λ) = gens(new_ring)
    F̂_regular = divexact(F̂_regular, ξ^valuation(F̂_regular, ξ))
    F¹_regular = obtain_evaluation(F̂_regular, α, 0)
    Pi¹_regular = Nemo.factor(F¹_regular)
    c = coeff(unit(Pi¹_regular), 1)
    @debug "An irreducible factorization at the first evaluation point is" Pi¹_regular

    # 3. Substitute ξ := 1
    F̂ = evaluate(F̂_regular, [ξ], [one(new_ring)])
    F¹ = evaluate(F¹_regular, [ξ], [one(new_ring)])
    squarefree_factors = collect(keys(Pi¹_regular.fac))
    powers = map(f -> Pi¹_regular.fac[f], squarefree_factors)
    Pi¹ = map(i -> evaluate(squarefree_factors[i]^powers[i], [ξ], [one(new_ring)]), 1:length(squarefree_factors))
    @debug "After substituting ξ := 1, the factors become" Pi¹ c

    # 4. Create buffers for storing the evaluations
    quad_ring, (_..., λ4) = PolynomialRing(ground_field, ["t", "u", "v", "λ"])
    collected_evaluations = [parent_ring_change(F¹, quad_ring)]
    collected_factor_evaluations = [map(f -> parent_ring_change(f, quad_ring), Pi¹)]

    # 5. Main evaluate-lift-interpolate loop
    m = 2
    Pi = nothing
    while true
        m0 = length(collected_evaluations) + 1
        @info "Computing the factorizations at indices $m0..$m"

        # Evaluate
        for i in m0:m
            Fⁱ = obtain_evaluation(F̂, α, i - 1)
            push!(collected_evaluations, parent_ring_change(Fⁱ, quad_ring))
        end

        # Lift
        for i in m0:m
            Fⁱ = collected_evaluations[i]
            Pi_quad = collected_factor_evaluations[i - 1]
            Pⁱi = hensel_lift_quad_variate(Fⁱ, c, Pi_quad)
            @assert prod(Pⁱi) == divexact(Fⁱ, c)
            @assert all(i -> evaluate(Pⁱi[i], [λ4], [zero(λ4)]) == evaluate(Pi_quad[i], [λ4], [one(λ4)]), 1:length(Pⁱi))
            push!(collected_factor_evaluations, Pⁱi)
        end

        # Try to interpolate
        success, Pi = try_to_interpolate_coefficients(new_ring, ring, c, collected_factor_evaluations, α)

        # Postprocessing passes
        Pi = map(P -> parent_ring_change(evaluate(P, vcat(one(ξ), xs, one(t), zero(u), zero(v), zero(λ))), ring), Pi)
        if any(iszero, Pi)
            m = 2*m
            continue
        end
        Pi = map(P -> divexact(P, leading_coefficient(P)), Pi)
        for i in 1:length(Pi)
            for x in gens(parent(Pi[i]))
                _, Pi[i] = remove(Pi[i], x)
            end
        end

        # Check if the interpolation was successful
        success = success && prod(Pi) * leading_coefficient(F) == F

        if success
            @info "Successful interpolation!"
            break
        end

        m = 2*m
    end

    # Postprocessing pass
    F = evaluate(F, gens(parent(F)) .* inv.(dilation))
    Pi = map(P -> evaluate(P, gens(parent(P)) .* inv.(dilation)), Pi)
    Pi = map(P -> divexact(P, leading_coefficient(P)), Pi)

    return leading_coefficient(F), append!(Pi, monomial_multiples)
end

# Same as above, but using the tri-variate substitution
function irreducible_factorization_2(F)
    # 0. Remove monomial multiples
    ring, ground_field = parent(F), base_ring(F)
    monomial_multiples = Vector{typeof(F)}()
    for x in gens(ring)
        v, F = remove(F, x)
        if !iszero(v)
            push!(monomial_multiples, x^v)
        end
    end
    isconstant(F) && return coeff(F, 1), monomial_multiples

    # 0. Apply a random dilation (or dilatation?..)
    dilation = map(_ -> rand(ground_field), 1:nvars(ring))
    F = evaluate(F, gens(ring) .* dilation)
    @info "Using random dilation $(gens(ring)) -> $(gens(ring) .* dilation)"

    # 1. Apply the transformation, together with a regularizing map x_i -> x_i ξ^w_n
    α = map(ground_field, Primes.nextprimes(2, nvars(ring)))
    β = map(_ -> rand(ground_field), 1:nvars(ring))
    γ = map(_ -> rand(ground_field), 1:nvars(ring))
    F̂_regular = homotopy_substitution_with_regularizing_transform_2(F, α, β, γ)
    @debug "After the substitution (factored for convenience):" Nemo.factor(F̂_regular)

    # 2. Compute a 3-variate factorization of the first specialization
    new_ring = parent(F̂_regular)
    (ξ, xs..., t, u, λ) = gens(new_ring)
    F̂_regular = divexact(F̂_regular, ξ^valuation(F̂_regular, ξ))
    F¹_regular = obtain_evaluation_2(F̂_regular, α, 0)
    Pi¹_regular = Nemo.factor(F¹_regular)
    c = coeff(unit(Pi¹_regular), 1)
    @debug "An irreducible factorization at the first evaluation point is" Pi¹_regular

    # 3. Substitute ξ := 1
    F̂ = evaluate(F̂_regular, [ξ], [one(new_ring)])
    F¹ = evaluate(F¹_regular, [ξ], [one(new_ring)])
    squarefree_factors = collect(keys(Pi¹_regular.fac))
    powers = map(f -> Pi¹_regular.fac[f], squarefree_factors)
    Pi¹ = map(i -> evaluate(squarefree_factors[i]^powers[i], [ξ], [one(new_ring)]), 1:length(squarefree_factors))
    @debug "After substituting ξ := 1, the factors become" Pi¹ c

    # 4. Create buffers for storing the evaluations
    tri_ring, (_..., λ3) = PolynomialRing(ground_field, ["t", "u", "λ"])
    collected_evaluations = [parent_ring_change(F¹, tri_ring)]
    collected_factor_evaluations = [map(f -> parent_ring_change(f, tri_ring), Pi¹)]

    # 5. Main evaluate-lift-interpolate loop
    m = 2
    Pi = nothing
    while true
        m0 = length(collected_evaluations) + 1
        @info "Computing the factorizations at indices $m0..$m"

        # Evaluate
        for i in m0:m
            Fⁱ = obtain_evaluation_2(F̂, α, i - 1)
            push!(collected_evaluations, parent_ring_change(Fⁱ, tri_ring))
        end

        # Lift
        for i in m0:m
            Fⁱ = collected_evaluations[i]
            Pi_tri = collected_factor_evaluations[i - 1]
            Pⁱi = hensel_lift_tri_variate(Fⁱ, c, Pi_tri)
            @assert prod(Pⁱi) == divexact(Fⁱ, c)
            @assert all(i -> evaluate(Pⁱi[i], [λ3], [zero(λ3)]) == evaluate(Pi_tri[i], [λ3], [one(λ3)]), 1:length(Pⁱi))
            push!(collected_factor_evaluations, Pⁱi)
        end

        # Try to interpolate
        success, Pi = try_to_interpolate_coefficients_2(new_ring, ring, c, collected_factor_evaluations, α)

        # Postprocessing passes
        Pi = map(P -> parent_ring_change(evaluate(P, vcat(one(ξ), xs, zero(t), zero(u), zero(λ))), ring), Pi)
        Pi = map(P -> evaluate(P, inv.(γ) .* gens(parent(P))), Pi)
        if any(iszero, Pi)
            m = 2*m
            continue
        end
        Pi = map(P -> divexact(P, leading_coefficient(P)), Pi)
        for i in 1:length(Pi)
            for x in gens(parent(Pi[i]))
                _, Pi[i] = remove(Pi[i], x)
            end
        end

        # Check if the interpolation was successful
        success = success && prod(Pi) * leading_coefficient(F) == F

        if success
            @info "Successful interpolation!"
            break
        end

        m = 2*m
    end

    # Postprocessing pass
    F = evaluate(F, gens(parent(F)) .* inv.(dilation))
    Pi = map(P -> evaluate(P, gens(parent(P)) .* inv.(dilation)), Pi)
    Pi = map(P -> divexact(P, leading_coefficient(P)), Pi)

    return leading_coefficient(F), append!(Pi, monomial_multiples)
end

############################################################
########################## Utils ###########################
############################################################

# Returns
#   F( 
#       ξ^w1 x1 ((1 - λ) + α1 λ) (t + β1 u + γ1 v) ),
#       ...,
#       ξ^wn xn ((1 - λ) + αn λ) (t + βn u + γn v) ),
#   )
function homotopy_substitution_with_regularizing_transform(F, α, β, γ)
    ring, ground_field = parent(F), base_ring(F)
    w = find_regularizing_weight(F)
    _, (ξ, xs..., t, u, v, λ) = PolynomialRing(
        ground_field, vcat("ξ", map(string, gens(ring)), "t", "u", "v", "λ")
    )
    F̂ = evaluate(F, (ξ .^ w) .* xs .* ((1 - λ) .+ α .* λ) .* (t .+ β .* u .+ γ .* v))
    @info """
    Input polynomial is in $(nvars(ring)) variables over $ground_field.
    Performing the substitution xᵢ -> ξʷⁱ xᵢ ((1 - λ) + αᵢλ) (t + βᵢ u + γᵢ v))
    Using the following evaluation points:
    α = $α
    β = $β
    γ = $γ
    And regularizing weight:
    w = $w"""
    return F̂
end

function homotopy_substitution_with_regularizing_transform_2(F, α, β, γ)
    ring, ground_field = parent(F), base_ring(F)
    w = find_regularizing_weight(F)
    _, (ξ, xs..., t, u, λ) = PolynomialRing(
        ground_field, vcat("ξ", map(string, gens(ring)), "t", "u", "λ")
    )
    F̂ = evaluate(F, (ξ .^ w) .* xs .* ((1 - λ) .+ α .* λ) .* (t .+ β .* u .+ γ))
    @info """
    Input polynomial is in $(nvars(ring)) variables over $ground_field.
    Performing the substitution xᵢ -> ξʷⁱ xᵢ ((1 - λ) + αᵢλ) (t + βᵢ u + γᵢ))
    Using the following evaluation points:
    α = $α
    β = $β
    γ = $γ
    And regularizing weight:
    w = $w"""
    return F̂
end

# Returns 
#   F̂(α1^i, ..., αn^i, u, v, t, λ)
function obtain_evaluation(F̂, α, i)
    ring = parent(F̂)
    @assert nvars(ring) == length(α) + 5
    (ξ, xs..., t, u, v, λ) = gens(ring)

    αⁱ = ring.(α .^ i)
    Fⁱ = evaluate(F̂, vcat(ξ, αⁱ, t, u, v, λ))

    return Fⁱ
end

function obtain_evaluation_2(F̂, α, i)
    ring = parent(F̂)
    @assert nvars(ring) == length(α) + 4
    (ξ, xs..., t, u, λ) = gens(ring)

    αⁱ = ring.(α .^ i)
    Fⁱ = evaluate(F̂, vcat(ξ, αⁱ, t, u, λ))

    return Fⁱ
end

function hensel_lift_quad_variate(Fⁱ, c, F_prev_factorization)
    ring = parent(Fⁱ)
    k = base_ring(ring)
    quadvariate_ring, (t4,u4,v4,λ4) = PolynomialRing(base_ring(ring), ["t","u","v","λ"])
    bivariate_ring, (t2, λ2) =  PolynomialRing(base_ring(ring), ["t","λ"])
    univariate_ring, (t1) = PolynomialRing(base_ring(ring), "t")
    interpolation_ring, (u2, v2) = PolynomialRing(base_ring(ring), ["u","v"])

    P_factors = F_prev_factorization

    F_quad = parent_ring_change(Fⁱ, quadvariate_ring)
    Ps_quad = map(P -> parent_ring_change(P, quadvariate_ring), P_factors)
    
    @debug "" F_quad c Ps_quad
    @assert c * prod(map(f -> evaluate(f, [λ4], [one(λ4)]), Ps_quad)) == mod(F_quad, λ4)

    # Apply a "random" shift
    subs = [t4, t4 + u4 + 12323, t4 + v4 + 1241, λ4]
    subs_rev = [t4, u4 - t4 - 12323, v4 - t4 - 1241, λ4]
    @assert evaluate(evaluate(F_quad, subs), subs_rev) == F_quad

    F_quad = evaluate(F_quad, subs)
    Ps_quad = map(P -> evaluate(P, subs), Ps_quad)

    T = (degree(F_quad, u4) + 1)*(degree(F_quad, v4) + 1)
    l = degree(F_quad, λ4) + 1
    interpolator = PrimesBenOrTiwari(interpolation_ring, T)
    point = startingpoint(interpolator)

    @debug "lifting with parameters" l T point
    
    points = [point .^i for i in 0:2T-1]
    factor_lifted_evaluations = []
    for (i, (p1,p2)) in enumerate(points)
        F_bi = evaluate(F_quad, [t2, p1, p2, λ2])
        Ps_uni = map(P -> evaluate(P, [t2, p1, p2, one(bivariate_ring)]), Ps_quad)
        
        @assert mod(F_bi, λ2) == c*prod(Ps_uni)
        @assert isone(content_in(F_bi, t2))
        
        Ps_lifted = hensel_lift_2(divexact(F_bi, c), Ps_uni, l, λ2)
        @assert divexact(F_bi, c) == prod(Ps_lifted)
        @assert all(i -> mod(Ps_lifted[i], λ2) == Ps_uni[i], 1:length(Ps_lifted))

        push!(factor_lifted_evaluations, Ps_lifted)
    end

    Ps_result = [zero(quadvariate_ring) for _ in 1:length(P_factors)]
    for i in 1:length(Ps_result)
        for (k, t) in enumerate(monomials(factor_lifted_evaluations[1][i]))
            evals = map(j -> coeff(factor_lifted_evaluations[j][i], t), 1:2T)
            interpolant = interpolate!(interpolator, points, evals)
            mm = evaluate(monomial(t, 1), [t4, λ4])
            cc = evaluate(interpolant, [u4, v4])
            Ps_result[i] += mm*cc
        end
    end

    Ps_result = map(P -> evaluate(P, subs_rev), Ps_result)
    @assert Fⁱ == c*prod(Ps_result)

    return Ps_result
end

function hensel_lift_tri_variate(Fⁱ, c, F_prev_factorization)
    ring = parent(Fⁱ)
    k = base_ring(ring)
    trivariate_ring, (t3,u3,λ3) = PolynomialRing(base_ring(ring), ["t","u","λ"])
    bivariate_ring, (t2, λ2) =  PolynomialRing(base_ring(ring), ["t","λ"])
    univariate_ring, (t1) = PolynomialRing(base_ring(ring), "t")
    interpolation_ring, (u1,) = PolynomialRing(base_ring(ring), ["u",])

    P_factors = F_prev_factorization

    F_tri = parent_ring_change(Fⁱ, trivariate_ring)
    Ps_tri = map(P -> parent_ring_change(P, trivariate_ring), P_factors)
    
    @debug "" F_tri c Ps_tri
    @assert c * prod(map(f -> evaluate(f, [λ3], [one(λ3)]), Ps_tri)) == mod(F_tri, λ3)

    # Apply a "random" shift
    subs = [t3, t3 + u3 + 12323, λ3]
    subs_rev = [t3, u3 - t3 - 12323, λ3]
    @assert evaluate(evaluate(F_tri, subs), subs_rev) == F_tri

    F_tri = evaluate(F_tri, subs)
    Ps_tri = map(P -> evaluate(P, subs), Ps_tri)

    T = degree(F_tri, u3) + 1
    l = degree(F_tri, λ3) + 1
    interpolator = PrimesBenOrTiwari(interpolation_ring, T)
    point = startingpoint(interpolator)

    @debug "lifting with parameters" l T point
    
    points = [point .^i for i in 0:2T-1]
    factor_lifted_evaluations = []
    for (i, (p1,)) in enumerate(points)
        F_bi = evaluate(F_tri, [t2, p1, λ2])
        Ps_uni = map(P -> evaluate(P, [t2, p1, one(bivariate_ring)]), Ps_tri)
        
        @assert mod(F_bi, λ2) == c*prod(Ps_uni)
        @assert isone(content_in(F_bi, t2))
        
        Ps_lifted = hensel_lift_2(divexact(F_bi, c), Ps_uni, l, λ2)
        @assert divexact(F_bi, c) == prod(Ps_lifted)
        @assert all(i -> mod(Ps_lifted[i], λ2) == Ps_uni[i], 1:length(Ps_lifted))

        push!(factor_lifted_evaluations, Ps_lifted)
    end

    Ps_result = [zero(trivariate_ring) for _ in 1:length(P_factors)]
    for i in 1:length(Ps_result)
        for (k, t) in enumerate(monomials(factor_lifted_evaluations[1][i]))
            evals = map(j -> coeff(factor_lifted_evaluations[j][i], t), 1:2T)
            interpolant = interpolate!(interpolator, points, evals)
            mm = evaluate(monomial(t, 1), [t3, λ3])
            cc = evaluate(interpolant, [u3,])
            Ps_result[i] += mm*cc
        end
    end

    Ps_result = map(P -> evaluate(P, subs_rev), Ps_result)
    @assert Fⁱ == c*prod(Ps_result)

    return Ps_result
end

function try_to_interpolate_coefficients(new_ring, interpolation_ring, c, collected_factor_evaluations, α)
    (ξ, xs..., t, u, v, λ) = gens(new_ring)
    
    interpolator = PrimesBenOrTiwari(interpolation_ring, length(collected_factor_evaluations) >> 1)
    @assert α == startingpoint(interpolator)
    points = [α .^ i for i in 0:length(collected_factor_evaluations)-1]

    @info "Trying to interpolate from $(length(collected_factor_evaluations)) evaluations at powers of α = $α"
    
    Ps_result = [zero(new_ring) for _ in 1:length(collected_factor_evaluations[1])]
    for i in 1:length(Ps_result)
        for (k, t) in enumerate(monomials(collected_factor_evaluations[1][i]))
            evals = c .* [coeff(collected_factor_evaluations[j][i], t) for j in 1:length(collected_factor_evaluations)]
            interpolant = interpolate!(interpolator, points, evals)
            mm = parent_ring_change(monomial(t, 1), new_ring)
            cc = evaluate(interpolant, xs)
            Ps_result[i] += mm*cc
        end
    end

    return true, Ps_result
end

function try_to_interpolate_coefficients_2(new_ring, interpolation_ring, c, collected_factor_evaluations, α)
    (ξ, xs..., t, u, λ) = gens(new_ring)
    
    interpolator = PrimesBenOrTiwari(interpolation_ring, length(collected_factor_evaluations) >> 1)
    @assert α == startingpoint(interpolator)
    points = [α .^ i for i in 0:length(collected_factor_evaluations)-1]

    @info "Trying to interpolate from $(length(collected_factor_evaluations)) evaluations at powers of α = $α"
    
    Ps_result = [zero(new_ring) for _ in 1:length(collected_factor_evaluations[1])]
    for i in 1:length(Ps_result)
        for (k, t) in enumerate(monomials(collected_factor_evaluations[1][i]))
            evals = c .* [coeff(collected_factor_evaluations[j][i], t) for j in 1:length(collected_factor_evaluations)]
            interpolant = interpolate!(interpolator, points, evals)
            mm = parent_ring_change(monomial(t, 1), new_ring)
            cc = evaluate(interpolant, xs)
            Ps_result[i] += mm*cc
        end
    end

    return true, Ps_result
end

############################################################
####################### Inline tests #######################
############################################################

using Nemo
using Test

R, (x1,x2,x3) = Nemo.GF(2^31-1)["x1","x2","x3"]

F = 4x1 + 4x2
c, Pi = irreducible_factorization(F)
@test c == 4*one(c) && Pi == [x1 + x2]

F = x1*x2
c, Pi = irreducible_factorization(F)
@test c == one(c) && Pi == [x1, x2]

F = (x1 + x2) * (x1*x2 + 1)
c, Pi = irreducible_factorization(F)
@test isone(c) && prod(Pi) * c == F

F = (x1 + x2) * (x1*x2 + 1) * (42x2 + 1)
c, Pi = irreducible_factorization(F)
@test c == 42*one(c) && prod(Pi) * c == F

F = (x1*x2 - x2*x3 + 1)^3*(x1*x2 + x2*x3)
c, Pi = irreducible_factorization(F)
@test prod(Pi) * c == F && length(Pi) == 3

F = (x1*x2 - x2*x3 + 1)*(x2 + x3)
c, Pi = irreducible_factorization(F)
@test prod(Pi) * c == F && length(Pi) == 2

F = (5x1*x2*x3 - 1) * (6x1 * x2 - 2) * (7x1 - 3)^2
c, Pi = irreducible_factorization(F)
@test prod(Pi) * c == F && length(Pi) == 3

F = (x1 + x2)*(x2 + x3)*(x1*x2 - x1*x3 - 2)*(x1*x3 - x2*x3 - 3)*(x1*x2*x3 - 99)
c, Pi = irreducible_factorization(F)
@test prod(Pi) * c == F && length(Pi) == 5

F = (x1 + x2)^2 + x3
c, Pi = irreducible_factorization(F)
@test prod(Pi) * c == F && length(Pi) == 1

F = x1^2 + x2^2 + x3^2
c, Pi = irreducible_factorization(F)
@test prod(Pi) * c == F && length(Pi) == 1

# Test tri-variate

F = 4x1 + 4x2
c, Pi = irreducible_factorization_2(F)
@test c == 4*one(c) && Pi == [x1 + x2]

F = x1*x2
c, Pi = irreducible_factorization_2(F)
@test c == one(c) && Pi == [x1, x2]

F = (x1 + x2) * (x1*x2 + 1)
c, Pi = irreducible_factorization_2(F)
@test isone(c) && prod(Pi) * c == F

F = (x1 + x2) * (x1*x2 + 1) * (42x2 + 1)
c, Pi = irreducible_factorization_2(F)
@test c == 42*one(c) && prod(Pi) * c == F

F = (x1*x2 - x2*x3 + 1)^3*(x1*x2 + x2*x3)
c, Pi = irreducible_factorization_2(F)
@test prod(Pi) * c == F && length(Pi) == 3

F = (x1*x2 - x2*x3 + 1)*(x2 + x3)
c, Pi = irreducible_factorization_2(F)
@test prod(Pi) * c == F && length(Pi) == 2

F = (5x1*x2*x3 - 1) * (6x1 * x2 - 2) * (7x1 - 3)^2
c, Pi = irreducible_factorization_2(F)
@test prod(Pi) * c == F && length(Pi) == 3

F = (x1 + x2)*(x2 + x3)*(x1*x2 - x1*x3 - 2)*(x1*x3 - x2*x3 - 3)*(x1*x2*x3 - 99)
c, Pi = irreducible_factorization_2(F)
@test prod(Pi) * c == F && length(Pi) == 5

F = (x1 + x2)^2 + x3
c, Pi = irreducible_factorization_2(F)
@test prod(Pi) * c == F && length(Pi) == 1

F = x1^2 + x2^2 + x3^2
c, Pi = irreducible_factorization_2(F)
@test prod(Pi) * c == F && length(Pi) == 1
