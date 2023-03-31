
mutable struct SharedStateJM{T,Ring,UnivRing}
    # univariate polynomial ring K[z]
    univring::UnivRing
    # multivariate polynomial ring K[x1, x2, ..., xn]
    ring::Ring
    # from 0 to n (outer loop of interpolation)
    varidx::Int
    # from 1 to 2t + 1 (inner loop of interpolation)
    tidx::Int
    tobeinitialized::Bool
    varchanged::Bool
    rootclash::Bool
    monomclash::Bool
    # the number of attempts to reinstate the interpolation
    # in case monomial evaluations clash (defaults to 0 or 1)
    attempts::Int
    # initial evaluation point,
    # α = (α1, α2, ..., αn), where αi is a generator of K
    α::Vector{T}
    # αn1 and αn2 are two generators of K
    αn1::T
    αn2::T
    base::Vector{T}
    # β is the current evaluation point
    β::Vector{T}
end

# the degree of the first entry in the geometric sequence
initialdegree() = 1

function SharedStateJM(ring::Ring) where {Ring}
    K = base_ring(ring)
    reinstateattempts = 1
    SharedStateJM(
        univariatize(AbstractAlgebra.PolyRing, ring),
        ring, 
        0, initialdegree(), 
        true, false, false, false, 
        reinstateattempts, 
        zeros(K, 0), zero(K), zero(K), zeros(K, 0), zeros(K, 0)
    )
end

mutable struct JavadiMonagan{T,Ring,UnivRing,UnivPoly} <: AbstractPolynomialInterpolator
    shared::SharedStateJM{T, Ring}
    bm::BerlekampMassey{UnivRing, T, UnivPoly}
    bmi::Vector{BerlekampMassey{UnivRing, T, UnivPoly}}
    L1::UnivPoly
    Li::Vector{UnivPoly}
    rs::Vector{T}
    rsi::Vector{Vector{T}}
    rstilde::Vector{Vector{T}}
    monoms::Vector{Vector{Int}}
    vs::Vector{T}
    actualt::Int
    boundd::Int
    isready::Bool
end

function JavadiMonagan(shared::SharedStateJM{T, Ring}, d::Integer) where {T, Ring}
    @assert d >= 0
    K = base_ring(shared.ring)
    n = nvars(shared.ring)
    bm = BerlekampMassey(shared.univring)
    bmi = map(c -> copy(bm), 1:n)
    JavadiMonagan(
        shared, 
        bm, 
        bmi, 
        zero(shared.univring), 
        Vector{elem_type(shared.univring)}(undef, n), 
        Vector{elem_type(K)}(undef, 0),
        Vector{Vector{elem_type(K)}}(undef, n), 
        Vector{Vector{elem_type(K)}}(undef, n),
        Vector{Vector{Int}}(undef, 0), 
        Vector{elem_type(K)}(undef, 0),
        0, d, 
        false
    )
end

mutable struct SimultaneousJavadiMonagan{T,Ring,UnivRing,UnivPoly,Poly}
    # a shared state of all interpolators
    shared::SharedStateJM{T, Ring}
    # a vector of interpolators
    jms::Vector{JavadiMonagan{T, Ring, UnivRing, UnivPoly}}
    # a vector of results as tuples (is interpolated?, polynomial)
    results::Vector{Tuple{Bool, Poly}}
end

function SimultaneousJavadiMonagan(ring::Ring, degrees::Vector{Int}) where {Ring}
    shared = SharedStateJM(ring)
    SimultaneousJavadiMonagan(shared, results)
end

function SimultaneousJavadiMonagan(shared::SharedStateJM, degrees::Vector{Int})
    jms = map(d -> JavadiMonagan(shared, d), degrees)
    results = map(_ -> (false, zero(shared.ring)), degrees)
    SimultaneousJavadiMonagan(shared, jms, results)
end

function allready(sjm::SimultaneousJavadiMonagan)
    all(x -> x[1], sjm.results)
end

function isready(sjm::SimultaneousJavadiMonagan, i::Integer)
    sjm.results[i][1]
end

function getresult(sjm::SimultaneousJavadiMonagan, i::Integer)
    sjm.results[i]
end

function nextpoint!(sjm::SimultaneousJavadiMonagan)
    if all(res -> res[1], sjm.results)
        @warn "All coefficients are already interpolated; doing nothing during this call."
        return sjm.shared.β
    end
    # update the shared state
    nextpoint!(sjm.shared)
end

function nextevaluation!(sjm::SimultaneousJavadiMonagan{T}, vs::Vector{T}) where {T}
    if all(res -> res[1], sjm.results)
        @warn "All coefficients are already interpolated; doing nothing during this call."
        return sjm.results
    end
    varchanged = true
    # advance the interpolators
    for i in 1:length(vs)
        status = nextevaluation!(sjm.jms[i], vs[i])
        # @info "" i status
        varchanged = varchanged && status.success
        if status.monomclash
            sjm.shared.monomclash = true
        end
    end
    if sjm.shared.monomclash
        if sjm.shared.attempts == 0
            @warn "No more attempts left"
            throw("No more attempts left")
        end
        return sjm.results
    end
    if varchanged
        sjm.shared.varidx += 1
        sjm.shared.varchanged = true
        sjm.shared.tidx = initialdegree()
    end
    for i in 1:length(vs)
        rootclash, idx = update!(sjm.jms[i], i, sjm.results)
        if rootclash
            for b in jm.bmi
                empty!(b)
            end
            sjm.shared.rootclash = true
            sjm.shared.varidx = 1
            sjm.shared.varchanged = false
            sjm.shared.tidx = initialdegree()
            return sjm.results
        end
    end
    sjm.results
end

function matchroots!(monoms, i, rs, rsbar, alphan, alphan1, dbound)
    t = length(rs)
    q = alphan1 // alphan[i]
    partial_degrees = zeros(Int, t)
    visitedi = fill(false, t)
    for j in 1:t
        d = 0
        success = false
        while d <= dbound && !success
            r = rs[j] * q^d
            for k in 1:t
                if !visitedi[k] && r == rsbar[k]
                    visitedi[k] = true
                    partial_degrees[j] = d
                    success = true
                    break
                end
            end
            d += 1
        end
        !success && return false
    end
    for j in 1:t
        monoms[j][i] = partial_degrees[j]
    end 
    return true
end

function matchroots!(monoms, i, rs, rsbar, rstilde, alphan, alphanp1, alphanp2, dbound)
    t = length(rs)
    q1 = alphanp1 // alphan[i]
    q2 = alphanp2 // alphan[i]
    # @info "" alphan alphanp1 alphanp2 q1 q2
    partial_degrees = zeros(Int, t)
    for j in 1:t
        d = 0
        success = false
        while d <= dbound && !success
            r1 = rs[j] * q1^d
            r2 = rs[j] * q2^d
            for k in 1:t
                if r1 == rsbar[k]
                    for kk in 1:t
                        if r2 == rstilde[kk]
                            partial_degrees[j] = d
                            success = true
                            break
                        end
                    end
                end
            end
            d += 1
        end
        @assert success
    end
    for j in 1:t
        monoms[j][i] = partial_degrees[j]
    end 
    return true
end

function nextpoint!(ssjm::SharedStateJM)
    R = ssjm.ring
    K = base_ring(R) 
    n = nvars(R)
    if ssjm.tobeinitialized
        # @info "tobeinitialized"
        alphas = distinct_generators(K, n + 1)
        alphan = alphas[1:n]
        alphanp1 = alphas[n + 1]
        alphanp2 = distinct_generator(K, alphan, alphanp1)
        # @info "" K n alphan alphanp1 alphanp2
        ssjm.α = alphan
        ssjm.αn1 = alphanp1
        ssjm.αn2 = alphanp2
        ssjm.base = alphan
        ssjm.varidx = 0
        ssjm.tidx = initialdegree()
        ssjm.tobeinitialized = false
    elseif ssjm.rootclash
        # @info "rootclash"
        ssjm.base = copy(ssjm.α)
        ssjm.base[ssjm.varidx] = ssjm.αn2
        # ssjm.rootclash = false
    elseif ssjm.monomclash
        # @info "monomclash"
        ssjm.monomclash = false
        ssjm.attempts -= 1
        ssjm.tobeinitialized = true
        return nextpoint!(ssjm)
    elseif ssjm.varchanged
        # @info "varchanged"
        ssjm.base = copy(ssjm.α)
        ssjm.base[ssjm.varidx] = ssjm.αn1
        ssjm.varchanged = false
    end
    ssjm.β = ssjm.base .^ ssjm.tidx
    # @info "" ssjm.β ssjm.base ssjm.tidx
    ssjm.tidx += 1
    ssjm.β
end

function nextevaluation!(jm::JavadiMonagan, v)
    varidx = jm.shared.varidx
    success = false
    monomclash = false
    if varidx == 0
        push!(jm.vs, v)
        success, L1 = next!(jm.bm, v)
        if success
            jm.L1 = Nemo.reverse(L1)
            jm.actualt = degree(jm.L1)
            # @info "success #0" jm.L1 jm.actualt
        end
    else
        success, Li = next!(jm.bmi[varidx], v)
        if success
            jm.Li[varidx] = reverse(Li)
            # @info "success #$varidx" jm.Li[varidx] jm.actualt
            if degree(Li) != jm.actualt
                @warn "Variable #$varidx: unlucky monomial clash in Javadi-Monagan."
                monomclash = true
                success = false
            end
        end
    end
    (success=success, monomclash=monomclash)
end

function update!(jm::JavadiMonagan, idx, results)
    R = jm.shared.ring
    K = base_ring(R)
    n = nvars(R)
    #
    !jm.shared.rootclash && !jm.shared.varchanged && return (false, 0)
    # 
    !jm.shared.rootclash && jm.shared.varidx <= n && return (false, 0)
    t = jm.actualt
    if !jm.shared.rootclash
        jm.rs = Nemo.roots(jm.L1)
        @assert length(jm.rs) == jm.actualt
        jm.monoms = map(_ -> zeros(Int, n), 1:t)
    end
    # 
    for i in 1:n
        if !jm.shared.rootclash
            jm.rsi[i] = Nemo.roots(jm.Li[i])
            @assert length(jm.rsi[i]) == t
            success = matchroots!(jm.monoms, i, jm.rs, jm.rsi[i], jm.shared.α,  jm.shared.αn1, jm.boundd)
            if !success
                @warn "Variable #$i: unlucky root clash in Javadi-Monagan."
                return (true, i)
            end
        else
            jm.rstilde[i] = Nemo.roots(jm.Li[i])
            success = matchroots!(jm.monoms, i, jm.rs, jm.rsi[i], jm.rstilde[i], jm.shared.α,  jm.shared.αn1, jm.shared.αn2, jm.boundd)
            @assert success
        end
    end
    Vspace = Nemo.MatrixSpace(K, t, t)
    A = Vspace(reshape([
        jm.rs[i] ^ (j)
        for i in 1:t
            for j in 1:t
    ], (t, t)))
    vi = jm.vs[1:t]
    ai = Nemo.inv(A) * vi
    results[idx] = true, R(ai, jm.monoms)
    (false, 0)
end

function interpolate!(jm::JavadiMonagan, blackbox)
    R = jm.ring
    K = base_ring(R) 
    n = nvars(R)
    attempts = 2
    L1 = zero(jm.bm.ring)
    @label Start
    while attempts > 0
        attempts -= 1
        empty!(jm.bm)
        bmi = [copy(jm.bm) for i in 1:n]
        alphas = distinct_generators(K, n + 1)
        alphan = alphas[1:n]
        alphanp1 = alphas[n + 1]
        alphanp2 = distinct_generator(K, alphan, alphanp1)
        # @info "" K n alphan alphanp1 alphanp2
        betas = Vector{Vector{elem_type(K)}}(undef, 0)
        vs = Vector{elem_type(K)}(undef, 0)
        success = false
        i = 1
        while !success
            beta = alphan .^ i
            v = blackbox(beta)
            i += 1
            push!(betas, beta)
            push!(vs, v)
            success, L1 = next!(jm.bm, v)
        end
        # @info "First BM done"
        L1 = Nemo.reverse(L1)
        Li = Vector{typeof(L1)}(undef, n)
        # the number of terms
        t = degree(L1)
        # # @info "" L1 t
        for i in 1:n
            for j in 1:2t + 1
                beta = copy(alphan)
                beta[i] = alphanp1
                beta = beta .^ j
                v = blackbox(beta)
                success, f = next!(bmi[i], v)
                Li[i] = reverse(f)
            end
            if degree(Li[i]) != t
                @warn "Javadi-Monagan failed due to monomial clash in variable $i"
                @goto Start
            end
        end
        rs = Nemo.roots(L1)
        # @info "BM done" rs
        @assert length(rs) == t
        rsi = Vector{typeof(rs)}(undef, n)
        monoms = map(_ -> zeros(Int, n), 1:t)
        for i in 1:n
            rsi[i] = Nemo.roots(Li[i])
            # @info "" i rsi[i]
            @assert length(rsi[i]) == t
            success, partial_degrees = matchroots(i, rs, rsi[i], alphan, alphanp1, jm.d)
            # @info "" success partial_degrees
            if !success
                @warn "Javadi-Monagan failed due to root clash"
                empty!(jm.bm)
                bmii = copy(jm.bm)
                Liitilde = zero(bmii.ring)
                for j in 1:2t + 1
                    beta = copy(alphan)
                    beta[i] = alphanp2
                    beta = beta .^ j
                    v = blackbox(beta)
                    success, f = next!(bmii, v)
                    Liitilde = reverse(f)
                end
                @assert success
                rsiitilde = Nemo.roots(Liitilde)
                # @info "" rsiitilde
                @assert length(rsiitilde) == t
                success, partial_degrees = matchroots(i, rs, rsi[i], rsiitilde, alphan, alphanp1, alphanp2, jm.d)
            end
            @assert success
            for j in 1:t
                monoms[j][i] = partial_degrees[j]
            end
        end
        Vspace = Nemo.MatrixSpace(K, t, t)
        A = Vspace(reshape([
            rs[i] ^ (j)
            for i in 1:t
                for j in 1:t
        ], (t, t)))
        vi = vs[1:t]
        ai = Nemo.inv(A) * vi
        return true, R(ai, monoms)
    end
    return false, zero(R)
end
