# Define the base field (actually an extension of F2)
R.<Z> = PolynomialRing(GF(2))
F2m.<z> = GF(2^127, modulus = Z^127 + Z^63 + 1)
# Define a quadratic extension on top)
F2x.<s> = F2m.extension(x^2 + x + 1)

def derive_trace_2m_fast_formula():
    print("z^i with nonzero trace in F_q")
    for i in range(127):
        elem = F2m(z^i)
        tr = elem.trace()
        if(tr != 0):
            print(i, tr)
    print()

derive_trace_2m_fast_formula()