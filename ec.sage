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

def curve_details(b):
    #Define a curve that E is a quadratic twist of
    q = ZZ(2^127)
    E = EllipticCurve(F2m, [1, 1, 0, 0, b])
    t = q + 1 - E.cardinality()
    #t=-t also gives solution
    n = (q-1)^2 + t^2
    
    r = max(n.factor())[0]
    S = Integers(r)
    #Formula from "A New Double Point Multiplication Method..."
    mu = S(q-1)*(S(t)^(-1))
    return n, r, t, mu

def curve_details_test():
    print("Based on current curve")
    n1, r1, t, mu = curve_details(b)

    assert n == n1
    assert r == r1
    assert mu^2 + 1 == 0

    print("Based on values in paper A New Double Point Multiplication Method")
    b1 = F2m.fetch_int(0x2BACF997126F185C3E67CB944EEB1168)
    n1, r1, t, mu = curve_details(b1)
    S = Integers(r1)
    assert n1 == 2*r1
    assert r1 == ZZ(14474011154664524427946373126085988481488994894707048965286167243381422079089)
    assert mu^2 + 1 == 0
    #assert S(mu) == S(8008021148421066531327005693257209127155969932631024546964258714431222018403)
    print("Curve details works!")

#As described by Karabina et al
def decomp(k, q, t):
    beta1 = (QQ(1-q)/QQ(t^2 + (q-1)^2))*QQ(k)
    beta2 = (QQ(t)/QQ(t^2 + (q-1)^2))*QQ(k)
    b1 = beta1.round()
    b2 = beta2.round()
    
    k1 = k - b1*(1-q) - b2*t
    k2 = - b1*t - b2*(q-1)
    return k1, k2

def decomp_test():
    n, r, t, mu = curve_details(b)
    k = ZZ(r-2)
    k1, k2 = decomp(k, ZZ(2^127), t)
    print(k1, k2)
    S = Integers(r)
    assert S(k1)+S(k2)*S(mu) == S(k)
    print("Decomp works!")

#Input must be odd
#Modified alg 6 from "Exponent Recoding and Regular Exponentiation Algorithms"
#by Joye et al
def regular_recode(k, w):
    l = ZZ(QQ(127/(w-1)).ceil()+1)
    v = ZZ(2^(w-1))
    acc = k
    k_reg = [0] * l
    
    for i in range(l-1):
        k_reg[i] = (n % (2*v))-v
        acc = (acc-k_reg[i])/v
        
    k_reg[l-1] = acc
    return k_reg

def regular_recode_test():
    k1 = ZZ(85070591730234615877113501116496779623)
    w = 4
    k1_reg = regular_recode(k1, w)
    assert(len(k1_reg) == 44)
    
    acc = ZZ(0)
    v = ZZ(1)
    for i in range(len(k1_reg)):
        acc = acc + k1_reg[i]*v
        v = v * 2^(w-1)
    assert acc == k1
    print("Regular recoding works!")

# We can see that the curve has a point (p, sqrt(b)) of small of small order 2
P = E(0, sqrt(b))
assert(2*P == 0*P)

# If this is a valid point in the curve, the following should be the point at infinity.
P = h * generator(E)
assert(r*P == 0*P)
Q = randrange(r) * P

curve_details_test()
decomp_test()
regular_recode_test()

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