
using Nemo

function construct_matrix(K, n)
    R, xs = PolynomialRing(K, ["x$i" for i in 1:n^2])
    S = MatrixSpace(R, n, n)  
    A = zero(S)
    for i in 1:n
        for j in 1:n
            A[i, j] = xs[(i-1)*n + j]
        end
    end
    A 
end

function toeplitz_matrix(K, n)
    R, xs = PolynomialRing(K, ["x$i" for i in 1:2n-1])
    S = MatrixSpace(R, n, n)   
    A = zero(S)
    for i in 1:n
        for j in i:n
            A[i, j] = xs[j-i+1]
        end
        for j in 1:i-1
            A[i, j] = xs[i - j + n]
        end
    end
    A 
end

function circulant_matrix(K, n)
    R, xs = PolynomialRing(K, ["x$i" for i in 1:n])
    S = MatrixSpace(R, n, n)   
    A = zero(S)
    for i in 1:n
        for j in 1:n
            A[i, j] = xs[abs((i - j) % n) + 1]
        end
    end
    A 
end


K = GF(2^31-1)
n = 4

A = circulant_matrix(K, 10);
d = det(A);
f = factor(d)