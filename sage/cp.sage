def compress(xP, lP):
    return xP + s*(lP[0].integer_representation() % 2)

def uncompress(cp):
    t = (cp[1].integer_representation() % 2) + F2m(1)
    xP = cp + t*s
    eq = b/xP^2 + xP^2 + a
    lP = half_trace_2x(eq)
    if (lP[0].integer_representation() % 2 != t):
        lP += 1
    return (xP, lP)
