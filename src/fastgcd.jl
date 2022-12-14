# Fast polynomial gcd algorithm.
# The notation and step numbering is taken from
# Algorithm 11.4, Modern Computer Algebra, by Gathen and Gerhard

# Given matrix A and vector x, returns the matrix-vector product Ax.
function matvec2by1(A, x)
    (A[1][1]*x[1] + A[1][2]*x[2], A[2][1]*x[1] + A[2][2]*x[2])
end

# Given two matrices A and B, returns the matrix product AB.
function matmul2by2(A, B)
    (A[1][1]*B[1][1] + A[1][2]*B[2][1], A[1][1]*B[1][2] + A[1][2]*B[2][2]),
    (A[2][1]*B[1][1] + A[2][2]*B[2][1], A[2][1]*B[1][2] + A[2][2]*B[2][2])
end

# Returns f ⥣ k, given as f quo x^(n - k) 
# (here, n = degree(f))
function ⥣(f, k)
    (k < 0) && return zero(f)
    if degree(f) < k
        Nemo.shift_left(f, k - degree(f))
    else
        Nemo.shift_right(f, degree(f) - k)
    end
end

# Fast polynomial gcd.
function fastgcd(r0, r1, k)
    n0, n1 = degree(r0), degree(r1)
    if iszero(r1) || k < n0 - n1
        # returns 0 and 
        # (1, 0)
        # (0, 1)
        return zero(r0), ((one(r0), zero(r0)), (zero(r0), one(r0)))
    end
    # first recursive call
    d = div(k, 2)
    jm1, R = fastgcd(r0 ⥣ 2d, r1 ⥣ (2d - (n0 - n1)), d)
    rjm1, rj = matvec2by1(R, (r0, r1))
    njm1, nj = degree(rjm1), degree(rj)
    if iszero(rj) || k < n0 - nj
        return jm1, R
    end
    qj = div(rjm1, rj)
    rjp1 = rjm1 - qj*rj
    rhojp1 = leading_coefficient(rjp1)
    # in Nemo, the leading_coefficient of 0 is 0;
    # we want it to be 1 here.
    iszero(rhojp1) && (rhojp1 = one(rhojp1))
    rjp1 = divexact(rjp1, rhojp1)
    njp1 = degree(rjp1)
    # second recursive call
    d⁺ = k - (n0 - nj)
    hmj, S = fastgcd(rj ⥣ 2d⁺, rjp1 ⥣ (2d⁺ - (nj - njp1)), d⁺)
    Qj = ((zero(r0), one(r0)), (inv(rhojp1), -qj*inv(rhojp1)))
    hmj + (jm1 + 1), matmul2by2(S, matmul2by2(Qj, R))
end

function standardize(g, f)
    !ismonic(g) && (g = divexact(g, leading_coefficient(g)))
    !ismonic(f) && (f = divexact(f, leading_coefficient(f)))
    if degree(g) < degree(f)
        g, f = f, g
    end
    if degree(g) == degree(f)
        g, f = f, g - f 
    end
    @assert degree(g) > degree(f)
    g, f
end

function fastgcd(g, f)
    (iszero(g) || iszero(f)) && (return one(g))
    g, f = standardize(g, f)
    _, R = fastgcd(g, f, degree(g))
    h = first(matvec2by1(R, (g, f)))
    divexact(h, leading_coefficient(h))
end


# given (polynomials) g and f (|f| >= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*g + s*f, |r| < k, where |r| is the maximal possible
function fastconstrainedEEA(g, f, k)
    @assert degree(g) > degree(f)
    # g, f = standardize(g, f)
    _, R = fastgcd(g, f, degree(g) - k - 1)
    ri, rj = matvec2by1(R, (g, f))
    if degree(ri) <= k
        t, s = R[1]
        return ri, t, s 
    else
        t, s = R[2]
        return rj, t, s
    end
end

function slowgcd(g, f)
    R = parent(g)  # = K[x]
    U = (one(R), zero(R), f)  # = (1, 0, f)
    V = (zero(R), one(R), g)  # = (0, 1, g)
    # in Nemo, degree(0) is -1
    while !iszero(V[3])
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end
    U[3]
end

# given (polynomials) g and f (|f| >= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*g + s*f, |r| < k, where |r| is the maximal possible
function constrainedEEA(g, f, k::Integer)
    @assert degree(f) >= degree(g)
    R = parent(g)  # = K[x]
    U = (one(R), zero(R), f)  # = (1, 0, f)
    V = (zero(R), one(R), g)  # = (0, 1, g)
    # in Nemo, degree(0) is -1
    while degree(V[3]) > k
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end
    (V[3], V[2], V[1])
end
