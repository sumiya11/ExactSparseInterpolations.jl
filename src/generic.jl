
const _bound = 2^20

generic_point(gf::Nemo.GaloisField) = rand(gf)
generic_point(gff::Nemo.GaloisFmpzField) = rand(gff)
generic_point(::Nemo.FlintIntegerRing) = ZZ(rand(1:_bound))
generic_point(::Nemo.FlintRationalField) = QQ(generic_point(ZZ))

function generic_point(ring)
    K = base_ring(ring)
    [generic_point(K) for _ in 1:nvars(ring)]
end

function geometric_point(ring)
    K = base_ring(ring)
    map(K, nextprimes(2, nvars(ring)))
end


univariatize(ring) = first(PolynomialRing(base_ring(ring), ["x"]))
true_univariatize(ring) = first(PolynomialRing(base_ring(ring), "x"))
