using Nemo

#####################################
# Generating random rest rational functions

# 
function random_rational_function(
        coeff_range, nvars,
        num_degree::Integer, den_degree::Integer,
        num_terms::Integer, den_terms::Integer; 
        ground_field=Nemo.QQ)
    R, _ = Nemo.PolynomialRing(ground_field, ["x$i" for i in 1:nvars])
    num = rand(R, num_terms:num_terms, 0:num_degree, coeff_range)
    den = rand(R, den_terms:den_terms, 0:den_degree, coeff_range)
    num // den
end

#####################################
# Test rational functions from
#   Interpolation of Dense and Sparse Rational
#   Functions and other Improvements in 
#   FireFly, 2021
# DOI: 10.1016/j.cpc.2021.107968

# Functions f1,f2,f3,f4 from Section 2.1.5 and Section 2.2
function FireFly_2021_testset1(; ground_field=Nemo.QQ)
    n1 = 20
    R1, zi1 = Nemo.PolynomialRing(ground_field, ["z$i" for i in 1:n1])
    f1 = sum(zi^20 for zi in zi1) // (sum(zi^20 for zi in zi1[1:10]) - sum(zi^20 for zi in zi1[11:20]))
    n2 = 5
    R2, (z1,z2,z3,z4,z5) = Nemo.PolynomialRing(ground_field, ["z$i" for i in 1:n2])
    f2 = (z1^100 + z2^200 + z3^300) // (z1*z2*z3*z4*z5 + (z1*z2*z3*z4*z5)^4)
    f3 = (1 + z1 + z2 + z3 + z4 + z5)^17 // (z4 - z2 + (z1*z2*z3*z4*z5)^10)
    f4 = ((1 + z1 + z2 + z3 + z4 + z5)^20 - 1) // (z4 - z2 + (z1*z2*z3*z4*z5)^10)
    [f1,f2,f3,f4]
end

# Functions f1,f2,f3,f4,f5 from Section 2.3
function FireFly_2021_testset2(;ground_field=Nemo.QQ)
    n1 = 4
    R1, (z1,z2,z3,z4) = Nemo.PolynomialRing(ground_field, ["z$i" for i in 1:n1])
    f1 = ((1 + z1 + z2 + z3 + z4 + z5)^20 - 1) // ((1 + z1 + z2 + z3 + z4 + z5)^20 - 1 + z^20)
    f2 = (z4 - 4)*(z4 - 3)*f1
    n2 = 20
    R2, zi2 = Nemo.PolynomialRing(ground_field, ["z$i" for i in 1:n2])
    f3 = sum(zi^20 for zi in zi1) // (sum(zi2[1]*zi2[2] + zi2[3]*zi2[4] + zi2[5]*zi2[6])^i for i in 1:5)
    f4 = f3 // (zi2[20]^35)
    f5 = f3 // (zi2[2]*zi2[20]^35 - 1)
    [f1,f2,f3,f4,f5]
end
