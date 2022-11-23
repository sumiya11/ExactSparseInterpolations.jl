
const _random_number_bound = 2^32-1

random_point(gf::Nemo.GaloisField) = rand(gf)
random_point(gff::Nemo.GaloisFmpzField) = rand(gff)
random_point(::Nemo.FlintIntegerRing) = ZZ(rand(1:_random_number_bound))
random_point(::Nemo.FlintRationalField) = QQ(random_point(ZZ))

function random_point(ring)
    K = base_ring(ring)
    map(_ -> random_point(K), 1:nvars(ring))
end

homogenize(ring; varname="x0") = first(PolynomialRing(base_ring(ring), append!([varname], map(string, AbstractAlgebra.symbols(ring)))))
dehomogenize(ring) = first(PolynomialRing(base_ring(ring), map(string, AbstractAlgebra.symbols(ring))[2:end]))

univariatize(::Type{Ring}, ring; varname="x") where {Ring<:AbstractAlgebra.MPolyRing} = first(PolynomialRing(base_ring(ring), [varname]))
univariatize(::Type{Ring}, ring; varname="x") where {Ring<:AbstractAlgebra.PolyRing} = first(PolynomialRing(base_ring(ring), varname))
