# Define the base field (actually an extension of F2)
R.<Z> = PolynomialRing(GF(2))
F2m.<z> = GF(2^127, modulus = Z^127 + Z^63 + 1)
# Define a quadratic extension on top)
F2x.<s> = F2m.extension(x^2 + x + 1)

def half_trace_2m(a):
    c = 0
    for i in range(0, 64):
        c = c + a^(2^(2*i))
    return c

def half_trace_2x(a):
    l1 = half_trace_2m(a[1])
    c0 = a[0] + l1^2
    l0 = half_trace_2m(c0 + c0.trace())
    return l0 + (l1 + c0.trace()) * s

def derive_trace_2m_fast_formula():
    print("z^i with nonzero trace in F_q")
    for i in range(127):
        elem = F2m(z^i)
        tr = elem.trace()
        if(tr != 0):
            print(i, tr)
    print()
