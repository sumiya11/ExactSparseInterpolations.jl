
# computes the product of all (z - x) for x in xs
# using the tree product algorithm.
# O(n log n) complexity, n = length(xs)
function producttree(z, xs)
    isempty(xs) && return one(z)
    length(xs) == 1 && return z - xs[1]
    n = length(xs)
    m = div(n, 2)
    producttree(z, view(xs, 1:m)) * producttree(z, view(xs, m+1:n))
end

treedepth(tree) = length(tree)
treebase(tree) = length(first(tree))

function buildproducttree(z::T, xs) where {T}
    n = length(xs)
    npow = nextpow(2, n)
    k = round(Int, log(2, npow)) + 1
    tree = Vector{Vector{T}}(undef, k)
    tree[1] = Vector{T}(undef, npow)
    for i in 1:n
        tree[1][i] = z - xs[i]
    end
    for i in n+1:npow
        tree[1][i] = one(z)
    end
    for i in 2:k
        nel = 2^(k-i)
        tree[i] = Vector{T}(undef, nel)
        for j in 1:nel
            tree[i][j] = tree[i-1][2*j-1] * tree[i-1][2*j]
        end
    end
    tree
end

function _remindertree!(f, rtree, ptree, depth, idx)
    iszero(f) && return rtree
    if depth == 0
        rtree[idx] = coeff(mod(f, ptree[1][idx]), 0)
        return rtree
    end
    l, r = 2*idx - 1, 2*idx
    r0 = mod(f, ptree[depth][l])
    r1 = mod(f, ptree[depth][r])
    _remindertree!(r0, rtree, ptree, depth - 1, l)
    _remindertree!(r1, rtree, ptree, depth - 1, r)
    rtree
end

function remindertree(f, ptree)
    K = base_ring(parent(f))
    rtree = zeros(K, treebase(ptree))
    _remindertree!(f, rtree, ptree, treedepth(ptree) - 1, 1)
end

function _lagrangetree(z, ys, ptree, depth, idx)
    R = parent(z)
    if depth == 0
        return R(ys[idx])
    end
    l, r = 2*idx - 1, 2*idx
    r0 = _lagrangetree(z, ys, ptree, depth - 1, l)
    r1 = _lagrangetree(z, ys, ptree, depth - 1, r)
    r0 * ptree[depth][r] + r1 * ptree[depth][l]
end

function lagrangetree(z, ys, ptree)
    _lagrangetree(z, ys, ptree, treedepth(ptree) - 1, 1) 
end

function fastpolyinterpolate(R, xs, ys)
    z = gen(R)
    ptree = buildproducttree(z, xs)
    m = ptree[end][1]
    dm = derivative(m)
    si = remindertree(dm, ptree)
    ysi = zeros(base_ring(R), length(ys))
    for i in 1:length(ys)
        ysi[i] = ys[i]*inv(si[i])
    end
    lagrangetree(z, ysi, ptree)
end
