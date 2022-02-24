# Define the base field (actually an extension of F2)
R.<Z> = PolynomialRing(GF(2))
F2m.<z> = GF(2^127, modulus = Z^127 + Z^63 + 1)
# Define a quadratic extension on top)
F2x.<s> = F2m.extension(x^2 + x + 1)

# Define curve coefficients, notice we are using a different b than the original paper.
a = F2x(s);
b = F2m(z^27 + 1)

E = EllipticCurve(F2x, [1, s, 0, 0, b])
# Sage does not compute the order of this curve, so I brought it over from MAGMA. It is 2 * r, for prime r.
n = 28948022309329048855892746252171976963485251808731388644510120425402211310058;
# Cofactor of the prime-order subgroup
h = 2
r = n//h

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

def generator(E):
    # Let's get a generator of the right order from MAGMA.
    x0 = F2m.fetch_int(0x5E0B72A98520F5A2D203CD2E4A5AE839)
    x1 = F2m.fetch_int(0x71B98581F8673A759639BBC43B8D797B)
    y0 = F2m.fetch_int(0x3ADACC9B694B43DB1D0CB95BEE9D4C31)
    y1 = F2m.fetch_int(0x3C8194E0263521C800C63FF2D65C6505)
    xP = x1*s + x0
    yP = y1*s + y0
    return E(xP,yP)

def to_lambda(P):
    return (P[0], P[0] + P[1]/P[0])

def compress(xP, lP):
    return (xP + s*lP[1].trace(), lP[0].trace())

def uncompress(cp):
    tr1 = 0
    xP, tr0 = cp
    if (xP[1].trace() != 1):
        # Fix bit of xP[1] and update the trace of lP[1]
        xP += s
        tr1 = xP[1].integer_representation() % 2
    t = b/xP^2 + xP^2 + a
    lP = half_trace_2x(t)
    if (lP[0].integer_representation() % 2 != tr0):
        lP += 1
    if (lP[1].integer_representation() % 2 != tr1):
        lP += s
    return (xP, lP)

# We can see that the curve has a point (p, sqrt(b)) of small of small order 2
P = E(0, sqrt(b))
assert(2*P == 0*P)

# If this is a valid point in the curve, the following should be the point at infinity.
P = generator(E)

for i in range(0, 100):
    # Pick a random point
    P = randrange(r) * P
    (xP, lP) = to_lambda(P)

    # Now let's check that the Weierstrass formulas to double a point are correct.

    # Now compute the formulas and verify, it should be true at the end.
    x2P = lP^2 + lP + a
    y2P = xP^2 + lP*x2P + x2P
    assert(E(x2P,y2P) == 2*P)

    # Check that lambda-coordinates satisfy equation
    assert((lP^2 + lP + a)*xP^2 == xP^4 + b)

    # Point compression method
    cp = compress(xP, lP)
    assert((xP, lP) == uncompress(cp))
