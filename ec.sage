load("cp.sage")
load("field.sage")

# Define curve coefficients, notice we are using a different b than the original paper.
one = F2m(1)
a = F2x(s);
b = F2m(z^27 + one)

E = EllipticCurve(F2x, [1, s, 0, 0, b])
# Sage does not compute the order of this curve, so I brought it over from MAGMA. It is 2 * r, for prime r.
n = 28948022309329048855892746252171976963485251808731388644510120425402211310058;
# Cofactor of the prime-order subgroup
h = 2
r = n//h

def generator(E):
    # Let's get a generator of the right order from MAGMA.
    x0 = F2m.fetch_int(0x5E0B72A98520F5A2D203CD2E4A5AE839)
    x1 = F2m.fetch_int(0x71B98581F8673A759639BBC43B8D797B)
    y0 = F2m.fetch_int(0x3ADACC9B694B43DB1D0CB95BEE9D4C31)
    y1 = F2m.fetch_int(0x3C8194E0263521C800C63FF2D65C6505)
    xP = x1*s + x0
    yP = y1*s + y0
    return E(xP,yP)

def to_lambda_aff(P):
    xP = P[0]
    lP = P[0] + P[1]/P[0]
    assert((lP^2 + lP + a)*xP^2 == xP^4 + b)
    return (xP, lP)

def to_lambda_prj(P):
    xP, lP = to_lambda_aff(P)
    return (xP, lP, one)

def from_lambda_aff(xP, lP):
    return E(xP, (lP + xP) * xP)

def from_lambda_prj(Xp, Lp, Zp):
    return E(Xp/Zp, (Lp/Zp + Xp/Zp) * Xp/Zp)

def double_weiss(xP, lP):
    # Now compute the formulas and verify, it should be true at the end.
    x2P = lP^2 + lP + a
    y2P = xP^2 + lP*x2P + x2P
    return (x2P, y2P)

def double_aff(xP, lP):
    x2P = lP^2 + lP + a
    l2P = xP^2 / x2P + lP^2 + a + 1
    return (x2P, l2P)

def double_prj(Xp, Lp, Zp):
    global mt, ma, mb, sq
    mt += 4
    ma += 1
    sq += 4
    T = Lp^2 + Lp * Zp + a * Zp^2
    X2 = T^2
    Z2 = T * Zp^2
    L2 = (Xp * Zp)^2 + X2 + T * (Lp * Zp) + Z2
    return (X2, L2, Z2)

def doubleb_prj(Xp, Lp, Zp):
    global mt, ma, mb, sq
    mt += 3
    ma += 1
    mb += 1
    sq += 4
    T = Lp^2 + Lp * Zp + a * Zp^2
    X2 = T^2
    Z2 = T * Zp^2
    L2 = (Lp + Xp)^2 * ((Lp + Xp)^2 + T + Zp^2) + (a^2 + b) * Zp^4 + X2 + (a + 1) * Z2
    return (X2, L2, Z2)

def add_prj(Xp, Lp, Zp, Xq, Lq, Zq):
    global mt, ma, mb, sq
    mt += 11
    sq += 2
    A = Lp * Zq + Lq * Zp
    B = (Xp * Zq + Xq * Zp)^2
    Xpq = A * (Xp * Zq) * (Xq * Zp) * A
    Lpq = (A * (Xq * Zp) + B)^2 + (A * B * Zq) * (Lp + Zp)
    Zpq = (A * B * Zq) * Zp
    return (Xpq, Lpq, Zpq)

def add_mix(Xp, Lp, Zp, xQ, lQ):
    global mt, ma, mb, sq
    mt += 8
    sq += 2
    A = Lp + lQ * Zp
    B = (Xp + xQ * Zp)^2
    Xpq = A * Xp * (xQ * Zp) * A
    Lpq = (A * (xQ * Zp) + B)^2 + (A * B) * (Lp + Zp)
    Zpq = (A * B) * Zp
    return (Xpq, Lpq, Zpq)

def smu_double_add(xP, lP, scalar):
    Xq = xP
    Lq = lP
    Zq = one
    for i in list(bin(k))[3:]:
        (Xq, Lq, Zq) = doubleb_prj(Xq, Lq, Zq)
        if int(i) != 0:
            (Xq, Lq, Zq) = add_prj(Xq, Lq, Zq, xP, lP, one)
    return (Xq, Lq, Zq)

def smu_double_always_add(xP, lP, scalar):
    (Xq, Lq, Zq) = (xP, lP, one)
    for i in list(bin(k))[3:]:
        (Xq, Lq, Zq) = doubleb_prj(Xq, Lq, Zq)
        (X3, L3, Z3) = add_prj(Xq, Lq, Zq, xP, lP, one)
        if int(i) != 0:
            (Xq, Lq, Zq) = (X3, L3, Z3)
    return (Xq, Lq, Zq)


# We can see that the curve has a point (p, sqrt(b)) of small of small order 2
P = E(0, sqrt(b))
assert(2*P == 0*P)

# If this is a valid point in the curve, the following should be the point at infinity.
P = h * generator(E)
assert(r*P == 0*P)
Q = randrange(r) * P

mt = ma = mb = sq = 0
for i in range(0, 10):
    k = randrange(r)

    # Pick a random point
    P = k * P
    Q = k * Q

    (xP, lP) = to_lambda_aff(P)
    (xQ, lQ) = to_lambda_aff(Q)

    # Test convertion to/from lambda coordinates
    assert(P == from_lambda_aff(xP, lP))
    assert(from_lambda_prj(xP, lP, one) == P)

    # Test point compression method
    cp = compress(xP, lP)
    assert((xP, lP) == uncompress(cp))

    # Test point doubling formulas
    assert(E(double_weiss(xP, lP)) == 2*P)

    (xP, lP) = double_aff(xP, lP)
    assert(from_lambda_aff(xP, lP) == 2*P)

    P = 2 * P
    (X2, L2, Z2) = double_prj(xP, lP, one)
    assert(from_lambda_prj(X2, L2, Z2) == 2*P)
    (X2, L2, Z2) = doubleb_prj(xP, lP, one)
    assert(from_lambda_prj(X2, L2, Z2) == 2*P)

    # Test point addition formulas
    (X3, L3, Z3) = add_mix(xP, lP, one, xQ, lQ)
    assert(from_lambda_prj(X3, L3, Z3) == P + Q)
    (X3, L3, Z3) = add_prj(xP, lP, one, xQ, lQ, one)
    assert(from_lambda_prj(X3, L3, Z3) == P + Q)

    # Test scalar multiplication algorithms
    (X3, L3, Z3) = smu_double_add(xP, lP, k)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_always_add(xP, lP, k)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)

# Benchmarks by operation counts
print("Operation counts as (muls, mul_a, mul_b, sqrs)")

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add(xP, lP, k)
print("Double-add       : ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_always_add(xP, lP, k)
print("Double-add-always: ", mt, ma, mb, sq)
