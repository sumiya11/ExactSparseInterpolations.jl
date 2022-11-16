#=
    Mahdi Jovadi and Monagan algorithm.

    Multivariate, Sparse, Las-Vegas randomized, 
    Blackbox

    Over integers modulo prime 
    
    4*n*T evaluations, where
        T is the upper bound on the number of terms,
        n is the number of variables
=#

mutable struct MahdiJovadiMonagan{T1,T2,T3,T4} <: AbstractInterpolator
    ring::T1
    D::Int
    t::Int
    T::Int
    bm::BerlekampMassey{T2, T3, T4}
end

function MahdiJovadiMonagan(ring, D::Integer, T::Integer; 
                t::Integer=-1,
                bm=BerlekampMassey(true_univariatize(ring)))
    MahdiJovadiMonagan(ring, D, t, T, bm)
end

function random_root_of_unity(Zp, Q)
    p = Int(characteristic(Zp))
    a = rand(Zp)
    while !root_of_unity(a)
        a = rand(Zp)
    end
    a
end

function root_of_unity(a::T) where {T}
    tmp = Set{T}(a^i for i in 0:Int(characteristic(parent(a))-1))
    length(tmp) == characteristic(parent(a)) - 1
end

function unique_random_roots_of_unity(n, Zp, Q)
    tmp = Set{elem_type(Zp)}()
    j = 0
    while length(tmp) < n
        j += 1
        j > 10^6 && error("too small field sorry")
        push!(tmp, random_root_of_unity(Zp, Q))
    end
    collect(tmp)
end

function evaluation_graph(r, r_k)
    return []
end

function has_unique_perfect_matching(G_k)
    return true
end

# over Zp for prime p
function Nemo.interpolate(mjm::MahdiJovadiMonagan{T1,T2,T3,T4}, blackbox) where {T1,T2,T3,T4}
    R = parent(mjm)
    Zp = base_ring(R)
    n = nvars(R)
    p = characteristic(Zp)
    T = mjm.T
    D = mjm.D
    @assert p > 0
    @assert true

    # label for restarting the algorithm
    @label FullRestart 
    empty!(mjm.bm)

    # generate n + 1 roots of unity in Zp
    Q = Primes.factor(Set, Int(p - 1))
    as = unique_random_roots_of_unity(n + 1, Zp, Q)
    alpha, alpha_n1 = as[1:n], as[n+1]
    @info "Roots of unity" alpha alpha_n1
    @assert all(root_of_unity, as)

    # generate 1 additional root, just in case
    alpha_n2 = alpha_n1
    gamma = alpha_n2
    while alpha_n2 in as
        gamma = random_root_of_unity(Zp, Q)
        alpha_n2 = gamma .* alpha_n2
    end
    @info "Gamma" gamma alpha_n2
    
    # prepare interpolation points
    beta = [alpha .^ i for i in 0:2T-1]
    
    # evaluate blackbox
    v = [blackbox(b) for b in beta]
    @info "Interpolation knots" beta v

    # find the generator of v0, v1, ...
    success, L = next!(mjm.bm, v[1])
    for i in 2:length(v)
        success, L = next!(mjm.bm, v[i])
    end
    @assert success
    # true number of nonzeros
    t = degree(L)
    @info "Recurrence generator" L t

    # compute the "integer" roots of generator
    r = Nemo.roots(L)
    @info "Generator roots" r

    # monomials to be intepolated
    monoms = [zeros(Int, n) for _ in 1:t]

    k = 1
    while k <= n
        @info "Variable $k.."
        # # determine degrees in variable xk:
        # # F := ∑ Ci*Mi^ei
        # # r[i] == Mi^ei
        # for i in 1:t
        #     di = 0
        #     Mi = alpha[k]
        #     while Mi^di != r[k]
        #         di += 1
        #     end
        #     monoms[t][k] = di
        # end

        # prepare new interpolation points
        # !!! take two additional points
        alpha_k = copy(alpha)
        alpha_k[k] = alpha_n1
        beta_k = [alpha_k .^ i for i in 0:2T+1]

        # evaluate blackbox
        v_k = [blackbox(b) for b in beta_k]
        @info "$k: " beta_k v_k

        # find the generator of v_k_0, v_k_1, ...
        empty!(mjm.bm)
        success, L_k = next!(mjm.bm, v_k[1])
        for i in 2:length(v_k)
            success, L_k = next!(mjm.bm, v_k[i])
        end
        @assert success
        @info "$k: interpolation knots" L_k

        # monomial evaluations are not distinct
        if degree(L_k) < t
            while alpha_n1 in as
                alpha_n1 = random_root_of_unity(Zp, Q)
            end
            @warn "Monomial collision. Selecting new generator"
            continue
        end

        # initial monomial evaluations collided
        if degree(L_k) > t
            @warn "Initial guess is not correct. Proceed with full restart"
            @goto FullRestart
        end

        # roots of new generator
        r_k = Nemo.roots(L_k)
        @info "$k: generator roots" r_k

        # construct the bipartite graph
        # between the roots of L and L_k
        G_k = evaluation_graph(r, r_k)

        if has_unique_perfect_matching(G_k)
            @info "Has the matching"
            # determine degrees in variable xk:
            # F := ∑ Ci*Mi^ei
            # r[i] == Mi^ei
            atoak = alpha_k[k] // alpha[k]
            for ri in r
                for di in 0:D
                    @info "" ri atoak di ri*atoak^di
                    if ri*atoak^di in r_k
                        monoms[t][k] = di
                        break
                    end
                end
            end
        else
            # ooops
        end

        k += 1
    end
    @info "Monoms" monoms

    # determine coefficients
    Vspace = Nemo.MatrixSpace(Zp, t, t)
    A = Vspace(reshape([
        r[i] ^ (j-1)
        for i in 1:t
            for j in 1:t
    ], (t, t)))

    vi = v[1:t]

    @info "" A vi
    
    ai = Nemo.inv(A) * vi

    R(ai, monoms)
end
