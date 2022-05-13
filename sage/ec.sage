from random import randrange


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

def add_psi(xP, lP):
    global mt, ma, mb, sq
    mt += 4
    sq += 1.5
    A = lP[1] + s
    B = xP[1]^2 # 1 subfield square
    Xpq = A * xP * A * (xP + xP[1]) # 4 subfield mults
    Zpq = (A * B) # 1 subfield mult
    Lpq = (A * (xP + xP[1]) + B)^2 + Zpq * (lP + 1) # 2 mults, 1 square
    return (Xpq, Lpq, Zpq)

def add_mix_mix(xP, lP, xQ, lQ):
    global mt, ma, mb, sq
    mt += 5
    sq += 2
    A = lP + lQ
    B = (xP + xQ)^2
    Xpq = A * xP * A * xQ
    Lpq = (A * xQ + B)^2 + (A * B) * (lP + 1)
    Zpq = (A * B)
    return (Xpq, Lpq, Zpq)

def add_sub_mix_mix(xP, lP, xQ, lQ):
    global mt, ma, mb, sq
    mt += 6
    sq += 4
    A = lP + lQ
    B = (xP + xQ)^2
    Xpq = A^2 * (xP * xQ)
    Zpq = (A * B)
    Lpq = (A * xQ + B)^2 + Zpq * (lP + 1)
    Xpmq = Xpq + xP*xQ
    Lpmq = Lpq + xQ^2 + B*(lP + 1)
    Zpmq = Zpq + B
    return (Xpq, Lpq, Zpq, Xpmq, Lpmq, Zpmq)

def double_add(Xq, Lq, Zq, xP, lP):
    global mt, ma, mb, sq
    mt += 10
    ma += 1
    sq += 6

    T = Lq^2 + Lq * Zq + a * Zq^2
    A = Xq^2 *Zq^2 + T*(Lq^2 + (a+1+lP)*Zq^2)
    B = (xP * Zq^2 + T)^2
    Xr = (xP * Zq^2) * A^2
    Zr = A * B *Zq^2
    Lr = T*(A+B)^2 + (lP+1)*Zr
    return (Xr, Lr, Zr)

def double_add_add(Xq, Lq, Zq, xP1, lP1, xP2, lP2):
    global mt, ma, mb, sq
    mt += 17
    ma += 1
    sq += 8

    T = Lq^2 + Lq*Zq + a*Zq^2
    U = Xq^2 * Zq^2 + T*(Lq^2 + (a + 1 + lP1)*Zq^2)
    F = xP1 * Zq^2
    G = (F+T)^2
    H = U^2 * F
    I = U * G * Zq^2
    J = (lP1 + lP2 + 1)*I + T*(U+G)^2

    K = xP2 * I
    L = (H + K)^2
    M = H * J
    Xr = J * K * M
    Zr = I * J * L
    Lr = (L + M)^2 + Zr*(lP2 + 1)
    return (Xr, Lr, Zr)

def neg_aff(xP, lP):
    return (xP, lP + 1)

def neg_proj(Xp, Lp, Zp):
    return (Xp, Lp + Zp, Zp)

def psi_aff(xP, lP):
    return (xP + xP[1], lP + lP[1] + s)

def curve_details(b):
    #Define a curve that E is a quadratic twist of
    q = ZZ(2^127)
    E = EllipticCurve(F2m, [1, 1, 0, 0, b])
    t = q + 1 - E.cardinality()
    #t=-t also gives solution
    n = (q-1)^2 + t^2

    r = n // 2 #r = max(n.factor())[0]
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
def decomp(k, t):
    q = ZZ(2^127)
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
    k1, k2 = decomp(k, t)
    print(k1, k2)
    S = Integers(r)
    assert S(k1)+S(k2)*S(mu) == S(k)
    print("Decomp works!")

def smu_double_add(xP, lP, scalar):
    (Xq, Lq, Zq) = (xP, lP, one)
    bits = list(reversed(Integer(scalar).bits()))
    for i in bits[1:]:
        (Xq, Lq, Zq) = doubleb_prj(Xq, Lq, Zq)
        if i:
            (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, xP, lP)
    return (Xq, Lq, Zq)

def smu_double_always_add(xP, lP, scalar):
    (Xq, Lq, Zq) = (xP, lP, one)
    bits = list(reversed(Integer(scalar).bits()))
    for i in bits[1:]:
        (Xq, Lq, Zq) = doubleb_prj(Xq, Lq, Zq)
        (X3, L3, Z3) = add_mix(Xq, Lq, Zq, xP, lP)
        if i: (Xq, Lq, Zq) = (X3, L3, Z3)
    return (Xq, Lq, Zq)

def smu_double_add_glv(xP, lP, scalar):
    n, r, t, mu = curve_details(b)
    k1, k2 = decomp(scalar, t)
    assert((k1 + k2*mu) % r == scalar)

    _xP, _lP = psi_aff(xP, lP)
    if (k1 < 0): lP += 1
    if (k2 < 0): _lP += 1

    bits1 = list(Integer(k1).bits())
    bits2 = list(Integer(k2).bits())

    l = max(len(bits1), len(bits2))
    if (len(bits1) > len(bits2)):
        (Xq, Lq, Zq) = (xP, lP, one)
    if (len(bits2) > len(bits1)):
        (Xq, Lq, Zq) = (_xP, _lP, one)
    if (len(bits1) == len(bits2)):
        (Xq, Lq, Zq) = add_mix(xP, lP, one, _xP, _lP)

    for i in range(l - 2, -1, -1):
        (Xq, Lq, Zq) = doubleb_prj(Xq, Lq, Zq)
        if (i < len(bits1) and bits1[i] != 0):
            (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, xP, lP)
        if (i < len(bits2) and bits2[i] != 0):
            (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, _xP, _lP)
    return (Xq, Lq, Zq)

def smu_double_add_glv_reg(xP, lP, scalar, w = 4):
    #b = F2m(z^49 + z^25 + 1)
    n, r, t, mu = curve_details(b)
    k1, k2 = decomp(scalar, t)
    assert((k1+k2*mu) % r == scalar)

    c1 = (k1 + 1) % 2
    c2 = (k2 + 1) % 2
    k1 = k1 + c1
    k2 = k2 + c2
    k1r = regular_recode(k1, w)
    k2r = regular_recode(k2, w)
    l = len(k1r)

    T = []
    (x2, l2) = double_aff(xP, lP)
    (Xacc, Lacc, Zacc) = (xP, lP, one)
    for i in range(2**(w-2)):
        T.append((Xacc / Zacc, Lacc / Zacc))
        (Xacc, Lacc, Zacc) = add_mix(Xacc, Lacc, Zacc, x2, l2)

    #Will convert table to affine coordinates using simultaneous inversion algorithm, so add costs here:
    global mt, ma, mb, sq
    mt += 5*(2**(w-2)-1) + 13

    _xP, _lP = psi_aff(xP, lP)
    (xP1, lP1, xP2, lP2) = (one, one, one, one)

    (xP1, lP1) = T[(abs(k1r[l-1])-1)/2]
    (xP2, lP2) = T[(abs(k2r[l-1])-1)/2]
    (xP2, lP2) = psi_aff(xP2, lP2)
    if k1r[l-1] < 0:
        (xP1, lP1) = neg_aff(xP1, lP1)
    if k2r[l-1] < 0:
        (xP2, lP2) = neg_aff(xP2, lP2)
    (Xq, Lq, Zq) = add_mix(xP1, lP1, one, xP2, lP2)

    for i in range(l - 2, -1, -1):
        for j in range(w-2):
            (Xq, Lq, Zq) = double_prj(Xq, Lq, Zq)
        (xP1, lP1) = T[(abs(k1r[i])-1)/2]
        (xP2, lP2) = T[(abs(k2r[i])-1)/2]
        (xP2, lP2) = psi_aff(xP2, lP2)

        if k1r[i] < 0:
            (xP1, lP1) = neg_aff(xP1, lP1)
        if k2r[i] < 0:
            (xP2, lP2) = neg_aff(xP2, lP2)
        (Xq, Lq, Zq) = double_add_add(Xq, Lq, Zq, xP1, lP1, xP2, lP2)

    if c1 == 1:
        (mxP, mlP) = neg_aff(xP, lP)
        (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, mxP, mlP)
    if c2 == 1:
        (mxP, mlP) = neg_aff(_xP, _lP)
        (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, mxP, mlP)

    return (Xq, Lq, Zq)

def smu_double_add_glv_reg_tab(xP, lP, scalar, w = 4):
    #b = F2m(z^49 + z^25 + 1)
    n, r, t, mu = curve_details(b)
    k1, k2 = decomp(scalar, t)
    assert((k1+k2*mu) % r == scalar)

    c1 = (k1 + 1) % 2
    c2 = (k2 + 1) % 2
    k1 = k1 + c1
    k2 = k2 + c2
    k1r = regular_recode(k1, w)
    k2r = regular_recode(k2, w)
    l = len(k1r)

    T1 = []
    T2 = []
    T  = []
    (x2, l2) = double_aff(xP, lP)
    (Xacc, Lacc, Zacc) = (xP, lP, one)
    for i in range(2**(w-2)):
        T1.append((Xacc / Zacc, Lacc / Zacc))
        (Xacc, Lacc, Zacc) = add_mix(Xacc, Lacc, Zacc, x2, l2)

    #Will convert table to affine coordinates using simultaneous inversion algorithm, so add costs here:
    global mt, ma, mb, sq
    mt += 5*(2**(w-2)-1) + 13

    T = [None] * (2**(w-2)) * (2**(w-2))
    for i in range(2**(w-2)):
        (xP1, lP1) = T1[i]
        (xP2, lP2) = psi_aff(xP1, lP1)
        (Xacc, Lacc, Zacc) = add_mix_mix(xP1, lP1, xP2, lP2)
        T[2**(w-2)*i+i] = (Xacc / Zacc, Lacc / Zacc)
        for j in range(i+1, 2**(w-2)):
            (xP2, lP2) = T1[j]
            (xP2, lP2) = psi_aff(xP2, lP2)
            (X3, L3, Z3, X4, L4, Z4) = add_sub_mix_mix(xP1, lP1, xP2, lP2)
            T[2**(w-2)*i+j] = (X3 / Z3, L3 / Z3)
            T[2**(w-2)*j+i] = psi_aff(X4 / Z4, L4 / Z4)

    #Will convert table to affine coordinates using simultaneous inversion algorithm, so add costs here:
    mt += 5*(2**(2*(w-2))-1) + 13

    _xP, _lP = psi_aff(xP, lP)
    (xP1, lP1, xP2, lP2) = (one, one, one, one)

    (xP1, lP1) = T1[(abs(k1r[l-1])-1)/2]
    (xP2, lP2) = T1[(abs(k2r[l-1])-1)/2]
    (xP2, lP2) = psi_aff(xP2, lP2)
    if k1r[l-1] < 0:
        (xP1, lP1) = neg_aff(xP1, lP1)
    if k2r[l-1] < 0:
        (xP2, lP2) = neg_aff(xP2, lP2)
    (Xq, Lq, Zq) = add_mix(xP1, lP1, one, xP2, lP2)

    for i in range(l - 2, -1, -1):
        for j in range(w-2):
            (Xq, Lq, Zq) = double_prj(Xq, Lq, Zq)

        if k1r[i] > 0 and k2r[i] > 0:
            k = 2**(w-2)*(abs(k1r[i])-1)/2 + (abs(k2r[i])-1)/2
            (xP1, lP1) = T[k]
        if k1r[i] < 0 and k2r[i] < 0:
            k = 2**(w-2)*(abs(k1r[i])-1)/2 + (abs(k2r[i])-1)/2
            (xP1, lP1) = T[k]
            (xP1, lP1) = neg_aff(xP1, lP1)
        if k1r[i] < 0 and k2r[i] > 0:
            k = 2**(w-2)*(abs(k2r[i])-1)/2 + (abs(k1r[i])-1)/2
            (xP1, lP1) = T[k]
            (xP1, lP1) = psi_aff(xP1, lP1)
        if k1r[i] > 0 and k2r[i] < 0:
            k = 2**(w-2)*(abs(k2r[i])-1)/2 + (abs(k1r[i])-1)/2
            (xP1, lP1) = T[k]
            (xP1, lP1) = psi_aff(xP1, lP1)
            (xP1, lP1) = neg_aff(xP1, lP1)
        (Xq, Lq, Zq) = double_add(Xq, Lq, Zq, xP1, lP1)

    if c1 == 1:
        (mxP, mlP) = neg_aff(xP, lP)
        (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, mxP, mlP)
    if c2 == 1:
        (mxP, mlP) = neg_aff(_xP, _lP)
        (Xq, Lq, Zq) = add_mix(Xq, Lq, Zq, mxP, mlP)

    return (Xq, Lq, Zq)

def smu_double_add_glv_reg_tab_precomp(xP, lP, w=4):
    T1 = []
    T2 = []
    T  = []
    (x2, l2) = double_aff(xP, lP)
    (Xacc, Lacc, Zacc) = (xP, lP, one)
    for i in range(2**(w-2)):
        T1.append((Xacc / Zacc, Lacc / Zacc))
        (Xacc, Lacc, Zacc) = add_mix(Xacc, Lacc, Zacc, x2, l2)

    #Will convert table to affine coordinates using simultaneous inversion algorithm, so add costs here:
    global mt, ma, mb, sq
    mt += 5*(2**(w-2)-1) + 13

    T = [None] * (2**(w-2)) * (2**(w-2))
    for i in range(2**(w-2)):
        (xP1, lP1) = T1[i]
        (xP2, lP2) = psi_aff(xP1, lP1)
        (Xacc, Lacc, Zacc) = add_mix_mix(xP1, lP1, xP2, lP2)
        T[2**(w-2)*i+i] = (Xacc / Zacc, Lacc / Zacc)
        for j in range(i+1, 2**(w-2)):
            (xP2, lP2) = T1[j]
            (xP2, lP2) = psi_aff(xP2, lP2)
            (X3, L3, Z3, X4, L4, Z4) = add_sub_mix_mix(xP1, lP1, xP2, lP2)
            T[2**(w-2)*i+j] = (X3 / Z3, L3 / Z3)
            T[2**(w-2)*j+i] = psi_aff(X4 / Z4, L4 / Z4)

    #Will convert table to affine coordinates using simultaneous inversion algorithm, so add costs here:
    mt += 5*(2**(2*(w-2))-1) + 13
    return T1, T


#Input must be odd
#Modified alg 6 from "Exponent Recoding and Regular Exponentiation Algorithms"
#by Joye et al
def regular_recode(k, w):
    l = ZZ(QQ(127/(w-1)).ceil()+1)
    v = ZZ(2^(w-1))
    acc = k
    k_reg = [0] * l

    for i in range(l-1):
        k_reg[i] = (acc % (2*v))-v
        acc = ZZ((acc-k_reg[i])/v)

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
R = randrange(r) * P

curve_details_test()
decomp_test()
regular_recode_test()

mt = ma = mb = sq = 0
_, _, _, mu = curve_details(b)

for i in range(0, 1):
    k = randrange(r)

    # Pick a random point
    P = k * P
    Q = k * Q
    R = k * R

    (xP, lP) = to_lambda_aff(P)
    (xQ, lQ) = to_lambda_aff(Q)
    (xR, lR) = to_lambda_aff(R)

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
    (X3, L3, Z3) = add_mix_mix(xP, lP, xQ, lQ)
    assert(from_lambda_prj(X3, L3, Z3) == P + Q)
    (X3, L3, Z3, X4, L4, Z4) = add_sub_mix_mix(xP, lP, xQ, lQ)
    assert(from_lambda_prj(X3, L3, Z3) == P + Q)
    assert(from_lambda_prj(X4, L4, Z4) == P - Q)
    (X3, L3, Z3) = add_psi(xP, lP)
    assert(from_lambda_prj(X3, L3, Z3) == P + int(mu)*P)

    #Test atomic formulas:
    (X3, L3, Z3) = double_add(xQ, lQ, one, xP, lP)
    assert(from_lambda_prj(X3, L3, Z3) == 2*Q + P)
    (X3, L3, Z3) = double_add_add(xQ, lQ, one, xP, lP, xR, lR)
    assert(from_lambda_prj(X3, L3, Z3) == 2*Q + P + R)

    # Test scalar multiplication algorithms
    (X3, L3, Z3) = smu_double_add(xP, lP, k)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_always_add(xP, lP, k)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv(xP, lP, k)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 4)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 5)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 6)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 7)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 4)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 5)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)
    (X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 6)
    assert(from_lambda_prj(X3, L3, Z3) == k*P)

    (xP, lP) = psi_aff(xP, lP)
    assert(from_lambda_aff(xP, lP) == int(mu)*P)

# Benchmarks by operation counts
print("Operation counts as (muls, mul_a, mul_b, sqrs)")

# mt = ma = mb = sq = 0
# (X3, L3, Z3) = smu_double_add(xP, lP, k)
# print("Double-add       : ", mt, ma, mb, sq)

# mt = ma = mb = sq = 0
# (X3, L3, Z3) = smu_double_always_add(xP, lP, k)
# print("Double-add-always: ", mt, ma, mb, sq)

# mt = ma = mb = sq = 0
# (X3, L3, Z3) = smu_double_add_glv(xP, lP, k)
# print("GLV-double-add: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 4)
print("GLV-reg-4-double-add: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 5)
print("GLV-reg-5-double-add: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg(xP, lP, k, 6)
print("GLV-reg-6-double-add: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 3)
print("GLV-reg-3-double-add-tab: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 4)
print("GLV-reg-4-double-add-tab: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 5)
print("GLV-reg-5-double-add-tab: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(X3, L3, Z3) = smu_double_add_glv_reg_tab(xP, lP, k, 6)
print("GLV-reg-6-double-add-tab: ", mt, ma, mb, sq)

mt = ma = mb = sq = 0
(T1, T) = smu_double_add_glv_reg_tab_precomp(xP, lP, 4)
print("GLV-reg-4-double-add-tab_precomp: ", mt, ma, mb, sq)
