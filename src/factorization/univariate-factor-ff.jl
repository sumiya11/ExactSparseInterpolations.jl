# Factoring polynomials over finite fields.
# In the notation of Modern Computer Algebra, by Gathen and Gerhard

# f^n mod g
function repeated_square_mod(f, n, g)
    iszero(n) && return one(f)
    isone(n) && return mod(f, g)
    @assert n > one(n)
    ffmodg = mod(f*f, g)
    if iseven(n)
        repeated_square_mod(ffmodg, div(n, 2), g)
    else
        mod(f * mod(repeated_square_mod(ffmodg, div(n-1, 2), g), g), g)
    end
end

function distinct_degree_factorization_ff(f)
    @assert issquarefree(f) && ismonic(f)
    R = parent(f)
    K, x = base_ring(R), gen(R)
    q = order(K)
    fi = f
    hi = x
    gi = elem_type(R)[]
    while !isone(fi)
        hi = repeated_square_mod(hi, q, f)
        g = gcd(hi - x, fi)
        fi = divexact(fi, g)
        push!(gi, g)
    end
    gi
end

# Returns a splitting polynomial of f for factors of degree d, 
# or failure
function equal_degree_splitting_ff(f, d)
    @assert issquarefree(f) && ismonic(f)
    R = parent(f)
    K, x = base_ring(R), gen(R)
    q = order(K)
    @assert isodd(q)
    n = degree(f)
    @assert iszero(mod(n, d))
    a = rand(R, 0:n - 1)
    g1 = gcd(a, f)
    !isone(g1) && return (true, g1)
    b = repeated_square_mod(a, div(q^d - 1, 2), f)
    g2 = gcd(b - one(b), f)
    !isone(g2) && !(g2 == f) && return (true, g2)
    (false, f)
end

function equal_degree_factorization_ff(f, d)
    @assert issquarefree(f) && ismonic(f)
    degree(f) == d && return [f]
    success, g = false, one(f)
    while !success
        success, g = equal_degree_splitting_ff(f, d)
    end
    ans1 = equal_degree_factorization_ff(g, d)
    ans2 = equal_degree_factorization_ff(divexact(f, g), d)
    append!(ans1, ans2)
end

function root_finding(f)
    R = parent(f)
    K, x = base_ring(R), gen(R)
    ans = Vector{elem_type(K)}()
    q = order(K)
    h = repeated_square_mod(x, q, f)
    g = gcd(h - x, f)
    isone(g) && return ans
    factors = equal_degree_factorization_ff(g, 1)
    map(f -> -coeff(f, 0), factors)
end

function factorization_ff(f::T) where {T}
    @assert issquarefree(f)
    R = parent(f)
    K, x = base_ring(R), gen(R)
    q = order(K)
    hi = x
    vi = divexact(f, leading_coefficient(f))
    U = Dict{T,Int}()
    i = 0
    while !isone(vi)
        i += 1
        hi = repeated_square_mod(hi, q, f)
        g = gcd(hi - x, vi)
        isone(g) && continue
        gs = equal_degree_factorization_ff(g, i)
        for gi in gs
            e = 0
            while true
                flag, vtmp = divides(vi, gi)
                !flag && break
                vi = vtmp
                e += 1
            end
            U[gi] = e
        end
    end
    U
end
