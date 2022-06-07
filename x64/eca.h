/* elliptic curve arithmetic */

void ec_print(__m128i px0, __m128i px1, __m128i pl0, __m128i pl1) {
    uint64_t u[2] = { 0 };
    uint64_t v[2] = { 0 };

    __m128i ONE = _mm_set_epi64x(0x1, 0x1);

    low_red_127_063_000(px0, px1, ONE);
	low_red_127_063_000(pl0, pl1, ONE);

    _mm_store_si128((__m128i *) u, px0);
    _mm_store_si128((__m128i *) v, px1);

    uint64_t x0[2] = { 0 };
    uint64_t x1[2] = { 0 };

    x0[0] = u[0];
    x0[1] = v[0];
    x1[0] = u[1];
    x1[1] = v[1];

    _mm_store_si128((__m128i *) u, pl0);
    _mm_store_si128((__m128i *) v, pl1);

    uint64_t l0[2] = { 0 };
    uint64_t l1[2] = { 0 };

    l0[0] = u[0];
    l0[1] = v[0];
    l1[0] = u[1];
    l1[1] = v[1];

    bn_print(x0, 2);
    bn_print(x1, 2);
    bn_print(l0, 2);
    bn_print(l1, 2);
}

void ec_print_prj(__m128i px0, __m128i px1, __m128i pl0, __m128i pl1, __m128i pz0, __m128i pz1) {
    uint64_t u[2] = { 0 };
    uint64_t v[2] = { 0 };

    ec_print(px0, px1, pl0, pl1);

    __m128i ONE = _mm_set_epi64x(0x1, 0x1);
    low_red_127_063_000(pz0, pz1, ONE);

    _mm_store_si128((__m128i *) u, pz0);
    _mm_store_si128((__m128i *) v, pz1);

    uint64_t z0[2] = { 0 };
    uint64_t z1[2] = { 0 };

    z0[0] = u[0];
    z0[1] = v[0];
    z1[0] = u[1];
    z1[1] = v[1];

    bn_print(z0, 2);
    bn_print(z1, 2);
}

/* [p > p] full doubling alternative */
void eca_dbl_mix(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1) {
    /* var */
    __m128i c0, c1, t0, t1;

    /* point doubling */
    /* T */
    low_sqr(&t0, &t1, pl0, pl1);
    t0 = _mm_xor_si128(t0, pl0); t1 = _mm_xor_si128(t1, pl1);
    t0 = _mm_xor_si128(t0, _mm_set_epi64x(0x1, 0));
    /* rz */
    *rz0 = t0; *rz1 = t1;
    /* rx */
    low_sqr(rx0, rx1, t0, t1);
    /* rl */
    low_sqr(rl0, rl1, px0, px1);
    low_mul(&c0, &c1, pl0, pl1, t0, t1);
    *rl0 = _mm_xor_si128(*rl0, c0); *rl1 = _mm_xor_si128(*rl1, c1);
    *rl0 = _mm_xor_si128(*rl0, *rz0); *rl1 = _mm_xor_si128(*rl1, *rz1);
    *rl0 = _mm_xor_si128(*rl0, *rx0); *rl1 = _mm_xor_si128(*rl1, *rx1);

    /* end */
    return;
}

/* [p > p] full doubling alternative */
void eca_dbl_ful(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i pz0, __m128i pz1) {
    /* var */
    __m128i b0, b1, c0, c1, t0, t1;
    __m128i x0, x1;

    /* point doubling */
    /* B */
    low_sqr(&b0, &b1, pz0, pz1);
    /* C */
    c0 = _mm_xor_si128(pl0, px0); c1 = _mm_xor_si128(pl1, px1);
    low_sqr(&c0, &c1, c0, c1);
    /* T */
    t0 = _mm_xor_si128(pl0, pz0); t1 = _mm_xor_si128(pl1, pz1);
    low_mul(&t0, &t1, t0, t1, pl0, pl1);
    low_mul_00u(&x0, &x1, b0, b1);
    t0 = _mm_xor_si128(t0, x0); t1 = _mm_xor_si128(t1, x1);

    /* rx */
    low_sqr(rx0, rx1, t0, t1);
    /* rz */
    low_mul(rz0, rz1, t0, t1, b0, b1);
    /* rl */
    *rl0 = _mm_xor_si128(c0, t0); *rl1 = _mm_xor_si128(c1, t1);
    *rl0 = _mm_xor_si128(*rl0, b0); *rl1 = _mm_xor_si128(*rl1, b1);
    low_mul(rl0, rl1, *rl0, *rl1, c0, c1);
    low_sqr(&x0, &x1, b0, b1);
    low_mul_27u(&x0, &x1, x0, x1);
    x0 = _mm_xor_si128(x0, *rx0); x1 = _mm_xor_si128(x1, *rx1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);
    low_mul_01u(&x0, &x1, *rz0, *rz1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);

    /* end */
    return;
}

/* [m > m] affine doubling */
void eca_dbl_aff(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1) {
    /* var */
    __m128i t0, t1, x0, x1;

    /* point doubling */
    /* rz */
    low_sqr(&t0, &t1, pl0, pl1);
    t0 = _mm_xor_si128(t0, pl0);
    t0 = _mm_xor_si128(t0, _mm_set_epi64x(0x1, 0x0));
    t1 = _mm_xor_si128(t1, pl1);
    /* rx */
    low_sqr(rx0, rx1, t0, t1);
    /* rl */
    low_sqr(&x0, &x1, px0, px1);
    *rl0 = _mm_xor_si128(x0, *rx0); *rl1 = _mm_xor_si128(x1, *rx1);
    low_mul(&x0, &x1, t0, t1, pl0, pl1);
    x0 = _mm_xor_si128(x0, t0); x1 = _mm_xor_si128(x1, t1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);
    low_inv_var(&t0, &t1, t0, t1);
	low_mul(rx0, rx1, *rx0, *rx1, t0, t1);
	low_mul(rl0, rl1, *rl0, *rl1, t0, t1);
	low_red_127_063_000(*rx0, *rx1, t0);
	low_red_127_063_000(*rl0, *rl1, t0);
}

/* [pp > p] full addition */
void eca_add_ful(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i pz0, __m128i pz1,
                 __m128i qx0, __m128i qx1,
                 __m128i ql0, __m128i ql1,
                 __m128i qz0, __m128i qz1) {
    /* var */
    __m128i a0, a1, b0, b1, c0, c1, d0, d1;
    __m128i s0, s1, t0, t1;
    __m128i x0, x1;

    /* point addition */
    /* S */
    low_mul(&s0, &s1, px0, px1, qz0, qz1);
    /* T */
    low_mul(&t0, &t1, qx0, qx1, pz0, pz1);

    /* A */
    low_mul(&a0, &a1, pl0, pl1, qz0, qz1);
    low_mul(&x0, &x1, ql0, ql1, pz0, pz1);
    a0 = _mm_xor_si128(a0, x0); a1 = _mm_xor_si128(a1, x1);
    /* B */
    b0 = _mm_xor_si128(s0, t0); b1 = _mm_xor_si128(s1, t1);
    low_sqr(&b0, &b1, b0, b1);
    /* C */
    low_mul(&c0, &c1, a0, a1, t0, t1);
    /* D */
    low_mul(&d0, &d1, a0, a1, b0, b1);
    low_mul(&d0, &d1, d0, d1, qz0, qz1);

    /* rx */
    low_mul(rx0, rx1, a0, a1, s0, s1);
    low_mul(rx0, rx1, *rx0, *rx1, c0, c1);
    /* rl */
    *rl0 = _mm_xor_si128(b0, c0); *rl1 = _mm_xor_si128(b1, c1);
    low_sqr(rl0, rl1, *rl0, *rl1);
    x0 = _mm_xor_si128(pl0, pz0); x1 = _mm_xor_si128(pl1, pz1);
    low_mul(&x0, &x1, x0, x1, d0, d1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);
    /* rz */
    low_mul(rz0, rz1, d0, d1, pz0, pz1);
}

/* [pm > p] mixed addition */
void eca_add_mix(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i pz0, __m128i pz1,
                 __m128i qx0, __m128i qx1,
                 __m128i ql0, __m128i ql1) {
    /* var */
    __m128i a0, a1, b0, b1, c0, c1, d0, d1;
    __m128i t0, t1, x0, x1;

    /* point addition */
    /* T */
    low_mul(&t0, &t1, qx0, qx1, pz0, pz1);

    /* A */
    low_mul(&a0, &a1, ql0, ql1, pz0, pz1);
    a0 = _mm_xor_si128(a0, pl0); a1 = _mm_xor_si128(a1, pl1);
    /* B */
    b0 = _mm_xor_si128(px0, t0); b1 = _mm_xor_si128(px1, t1);
    low_sqr(&b0, &b1, b0, b1);
    /* C */
    low_mul(&c0, &c1, a0, a1, t0, t1);
    /* D */
    low_mul(&d0, &d1, a0, a1, b0, b1);

    /* rx */
    low_mul(rx0, rx1, a0, a1, c0, c1);
    low_mul(rx0, rx1, *rx0, *rx1, px0, px1);
    /* rl */
    *rl0 = _mm_xor_si128(c0, b0); *rl1 = _mm_xor_si128(c1, b1);
    low_sqr(rl0, rl1, *rl0, *rl1);
    x0 = _mm_xor_si128(pl0, pz0); x1 = _mm_xor_si128(pl1, pz1);
    low_mul(&x0, &x1, x0, x1, d0, d1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);
    /* rz */
    low_mul(rz0, rz1, d0, d1, pz0, pz1);

    /* end */
    return;
}

/* [pm > p] mixed addition */
void eca_add_mix_complete(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i pz0, __m128i pz1,
                 __m128i qx0, __m128i qx1,
                 __m128i ql0, __m128i ql1) {
    /* var */
    __m128i a0, a1, b0, b1, c0, c1, d0, d1, e0, e1, msk0, msk1, msk2;
    __m128i sx0, sx1, sl0, sl1, sz0, sz1;
    __m128i ZERO = _mm_setzero_si128(), ONE = _mm_set_epi64x(0x0, 0x1);

    /* (Xq2, Lq2, Zq2) = double_mix(xQ, lQ) */
    eca_dbl_ful(&a0, &a1, &b0, &b1, &c0, &c1, qx0, qx1, ql0, ql1, ONE, ZERO);
    /* (Xpq, Lpq, Zpq) = add_mix(Xp, Lp, Zp, xQ, lQ) */
    eca_add_mix(&sx0, &sx1, &sl0, &sl1, &sz0, &sz1, px0, px1, pl0, pl1, pz0, pz1, qx0, qx1, ql0, ql1);
    /* Xq = xQ * Zp, Lq = lQ * Zp */
    low_mul(&d0, &d1, qx0, qx1, pz0, pz1);
    low_mul(&e0, &e1, ql0, ql1, pz0, pz1);

    /* if Zp == 0, return (xQ, lQ, 1) */
    msk0 = _mm_and_si128(_mm_cmpeq_epi64(pz0, ZERO), _mm_cmpeq_epi64(pz1, ZERO));
    /* if Xp == Xq and Lp == Lq, return (Xq2, Lq2, Zq2) */
    msk1 = _mm_and_si128(_mm_and_si128(_mm_cmpeq_epi64(d0, px0), _mm_cmpeq_epi64(d1, px1)),
                         _mm_and_si128(_mm_cmpeq_epi64(e0, pl0), _mm_cmpeq_epi64(e1, pl1)));
    msk2 = _mm_or_si128(msk0, msk1);

    *rx0 = *rx1 = *rl0 = *rl1 = *rz0 = *rz1 = ZERO;
    *rx0 = _mm_xor_si128(*rx0, _mm_and_si128(qx0, msk0));
    *rx1 = _mm_xor_si128(*rx1, _mm_and_si128(qx1, msk0));
    *rl0 = _mm_xor_si128(*rl0, _mm_and_si128(ql0, msk0));
    *rl1 = _mm_xor_si128(*rl1, _mm_and_si128(ql1, msk0));
    *rz0 = _mm_xor_si128(*rz0, _mm_and_si128(ONE, msk0));
    *rz1 = _mm_xor_si128(*rz1, _mm_and_si128(ZERO, msk0));

    *rx0 = _mm_xor_si128(*rx0, _mm_and_si128(a0, msk1));
    *rx1 = _mm_xor_si128(*rx1, _mm_and_si128(a1, msk1));
    *rl0 = _mm_xor_si128(*rl0, _mm_and_si128(b0, msk1));
    *rl1 = _mm_xor_si128(*rl1, _mm_and_si128(b1, msk1));
    *rz0 = _mm_xor_si128(*rz0, _mm_and_si128(c0, msk1));
    *rz1 = _mm_xor_si128(*rz1, _mm_and_si128(c1, msk1));

    /* return (Xpq, Lpq, Zpq) */
    *rx0 = _mm_xor_si128(*rx0, _mm_andnot_si128(msk2, sx0));
    *rx1 = _mm_xor_si128(*rx1, _mm_andnot_si128(msk2, sx1));
    *rl0 = _mm_xor_si128(*rl0, _mm_andnot_si128(msk2, sl0));
    *rl1 = _mm_xor_si128(*rl1, _mm_andnot_si128(msk2, sl1));
    *rz0 = _mm_xor_si128(*rz0, _mm_andnot_si128(msk2, sz0));
    *rz1 = _mm_xor_si128(*rz1, _mm_andnot_si128(msk2, sz1));

    /* end */
    return;
}

/* [mm> p] mixed-mixed addition */
void eca_add_mma(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i qx0, __m128i qx1,
                 __m128i ql0, __m128i ql1) {
    /* var */
    __m128i a0, a1, b0, b1, c0, c1;
    __m128i x0, x1;

    /* point addition */
    /* A */
    a0 = _mm_xor_si128(pl0, ql0); a1 = _mm_xor_si128(pl1, ql1);
    /* B */
    b0 = _mm_xor_si128(px0, qx0); b1 = _mm_xor_si128(px1, qx1);
    low_sqr(&b0, &b1, b0, b1);
    /* C */
    low_mul(&c0, &c1, a0, a1, qx0, qx1);

    /* rz */
    low_mul(rz0, rz1, a0, a1, b0, b1);
    /* rx */
    low_mul(rx0, rx1, c0, c1, a0, a1);
    low_mul(rx0, rx1, *rx0, *rx1, px0, px1);
    /* rl */
    *rl0 = _mm_xor_si128(c0, b0); *rl1 = _mm_xor_si128(c1, b1);
    low_sqr(rl0, rl1, *rl0, *rl1);
    low_mul(&x0, &x1, *rz0, *rz1, pl0, pl1);
    x0 = _mm_xor_si128(x0, *rz0); x1 = _mm_xor_si128(x1, *rz1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);

    /* end */
    return;
}

/* [mm > p] doubling-and-mixed-addition */
void eca_dbl_mma(__m128i *rx0, __m128i *rx1,
                     __m128i *rl0, __m128i *rl1,
                     __m128i *rz0, __m128i *rz1,
                     __m128i px0, __m128i px1,
                     __m128i pl0, __m128i pl1,
                     __m128i qx0, __m128i qx1,
                     __m128i ql0, __m128i ql1) {
    /* var */
    __m128i d0, d1, e0, e1;
    __m128i t0, t1, a0, a1, b0, b1;
    __m128i x0, x1;
    __m128i ONE;

    /* init */
    ONE = _mm_set_epi64x(0x0, 0x1);

    /* D */
    low_sqr(&d0, &d1, pl0, pl1);
    /* E */
    e0 = _mm_xor_si128(ql0, ONE); e1 = ql1;

    /* T */
    ONE = _mm_slli_si128(ONE, 8);
    t0 = _mm_xor_si128(d0, pl0); t1 = _mm_xor_si128(d1, pl1);
    t0 = _mm_xor_si128(t0, ONE);
    /* A */
    low_sqr(&a0, &a1, px0, px1);
    x0 = _mm_xor_si128(d0, e0); x1 = _mm_xor_si128(d1, e1);
    x0 = _mm_xor_si128(x0, ONE);
    low_mul(&x0, &x1, x0, x1, t0, t1);
    a0 = _mm_xor_si128(a0, x0); a1 = _mm_xor_si128(a1, x1);
    /* B */
    b0 = _mm_xor_si128(qx0, t0); b1 = _mm_xor_si128(qx1, t1);
    low_sqr(&b0, &b1, b0, b1);

    /* rx */
    low_sqr(rx0, rx1, a0, a1);
    low_mul(rx0, rx1, *rx0, *rx1, qx0, qx1);
    /* rz */
    low_mul(rz0, rz1, a0, a1, b0, b1);
    /* rl */
    *rl0 = _mm_xor_si128(a0, b0); *rl1 = _mm_xor_si128(a1, b1);
    low_sqr(rl0, rl1, *rl0, *rl1);
    low_mul(rl0, rl1, *rl0, *rl1, t0, t1);
    low_mul(&x0, &x1, e0, e1, *rz0, *rz1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);

    /* end */
    return;
}

/* [pm > p] doubling-and-addition */
void eca_dbl_add(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i pz0, __m128i pz1,
                 __m128i qx0, __m128i qx1,
                 __m128i ql0, __m128i ql1) {
    /* var */
    __m128i d0, d1, e0, e1, f0, f1, s0, s1;
    __m128i t0, t1, a0, a1, b0, b1, c0, c1;
    __m128i x0, x1;
    __m128i ONE;

    /* init */
    ONE = _mm_set_epi64x(0x0, 0x1);

    /* point doubling-and-addition */
    /* D */
    low_sqr(&d0, &d1, pl0, pl1);
    /* E */
    low_sqr(&e0, &e1, pz0, pz1);
    /* F */
    f0 = _mm_xor_si128(ql0, ONE); f1 = ql1;
    /* S */
    low_mul(&s0, &s1, qx0, qx1, e0, e1);

    /* T */
    low_mul(&t0, &t1, pl0, pl1, pz0, pz1);
    t0 = _mm_xor_si128(t0, d0); t1 = _mm_xor_si128(t1, d1);
    low_mul_00u(&x0, &x1, e0, e1);
    t0 = _mm_xor_si128(t0, x0); t1 = _mm_xor_si128(t1, x1);
    /* A */
    low_sqr(&a0, &a1, px0, px1);
    low_mul(&a0, &a1, a0, a1, e0, e1);
    x0 = _mm_xor_si128(f0, _mm_slli_si128(ONE, 8)); x1 = f1;
    low_mul(&x0, &x1, x0, x1, e0, e1);
    x0 = _mm_xor_si128(x0, d0); x1 = _mm_xor_si128(x1, d1);
    low_mul(&x0, &x1, x0, x1, t0, t1);
    a0 = _mm_xor_si128(a0, x0); a1 = _mm_xor_si128(a1, x1);
    /* B */
    b0 = _mm_xor_si128(s0, t0); b1 = _mm_xor_si128(s1, t1);
    low_sqr(&b0, &b1, b0, b1);

    /* rx */
    low_sqr(rx0, rx1, a0, a1);
    low_mul(rx0, rx1, *rx0, *rx1, s0, s1);
    /* rz */
    low_mul(rz0, rz1, a0, a1, b0, b1);
    low_mul(rz0, rz1, *rz0, *rz1, e0, e1);
    /* rl */
    *rl0 = _mm_xor_si128(a0, b0); *rl1 = _mm_xor_si128(a1, b1);
    low_sqr(rl0, rl1, *rl0, *rl1);
    low_mul(rl0, rl1, *rl0, *rl1, t0, t1);
    low_mul(&x0, &x1, f0, f1, *rz0, *rz1);
    *rl0 = _mm_xor_si128(*rl0, x0); *rl1 = _mm_xor_si128(*rl1, x1);

    /* end */
    return;
}

/* [pm > p] doubling-and-addition-subtraction */
void eca_dbl_add_sub(__m128i *sx0, __m128i *sx1,
                    __m128i *sl0, __m128i *sl1,
                    __m128i *sz0, __m128i *sz1,
                    __m128i *rx0, __m128i *rx1,
                    __m128i *rl0, __m128i *rl1,
                    __m128i *rz0, __m128i *rz1,
                    __m128i px0, __m128i px1,
                    __m128i pl0, __m128i pl1,
                    __m128i pz0, __m128i pz1,
                    __m128i qx0, __m128i qx1,
                    __m128i ql0, __m128i ql1) {
    /* var */
    __m128i a0, a1, b0, b1, c0, c1, t0, t1, x0, x1;
    __m128i ONE;

    /* point doubling */
    /* B */
    low_sqr(&b0, &b1, pz0, pz1);
    /* C */
    c0 = _mm_xor_si128(pl0, px0); c1 = _mm_xor_si128(pl1, px1);
    low_sqr(&c0, &c1, c0, c1);
    /* T */
    t0 = _mm_xor_si128(pl0, pz0); t1 = _mm_xor_si128(pl1, pz1);
    low_mul(&t0, &t1, t0, t1, pl0, pl1);
    low_mul_00u(&x0, &x1, b0, b1);
    t0 = _mm_xor_si128(t0, x0); t1 = _mm_xor_si128(t1, x1);

    /* sx */
    low_sqr(sx0, sx1, t0, t1);
    /* sz */
    low_mul(sz0, sz1, t0, t1, b0, b1);
    /* sl */
    *sl0 = _mm_xor_si128(c0, t0); *sl1 = _mm_xor_si128(c1, t1);
    *sl0 = _mm_xor_si128(*sl0, b0); *sl1 = _mm_xor_si128(*sl1, b1);
    low_mul(sl0, sl1, *sl0, *sl1, c0, c1);
    low_sqr(&x0, &x1, b0, b1);
    low_mul_27u(&x0, &x1, x0, x1);
    x0 = _mm_xor_si128(x0, *sx0); x1 = _mm_xor_si128(x1, *sx1);
    *sl0 = _mm_xor_si128(*sl0, x0); *sl1 = _mm_xor_si128(*sl1, x1);
    low_mul_01u(&x0, &x1, *sz0, *sz1);
    *sl0 = _mm_xor_si128(*sl0, x0); *sl1 = _mm_xor_si128(*sl1, x1);

    /* A = lP * Z2 + L2, B' = xP*Z2. */
    low_mul(&a0, &a1, ql0, ql1, *sz0, *sz1);
    low_mul(&b0, &b1, qx0, qx1, *sz0, *sz1);
    a0 = _mm_xor_si128(a0, *sl0);
    a1 = _mm_xor_si128(a1, *sl1);
    /* C = (B' * X2) */
    low_mul(&c0, &c1, b0, b1, *sx0, *sx1);
    /* B = (B' + X2)^2 */
    b0 = _mm_xor_si128(b0, *sx0);
    b1 = _mm_xor_si128(b1, *sx1);
    low_sqr(&b0, &b1, b0, b1);

    ONE = _mm_set_epi64x(0x0, 0x1);
    /* Xpq = A^2 * C. */
    low_sqr(&t0, &t1, a0, a1);
    low_mul(rx0, rx1, t0, t1, c0, c1);
    /* Zpq = A * B * Z2 */
    low_mul(rz0, rz1, *sz0, *sz1, a0, a1);
    low_mul(rz0, rz1, *rz0, *rz1, b0, b1);
    /* Lpq = ((A*X2) + B)^2 + Zpq*(lP + 1) */
    low_mul(&t0, &t1, a0, a1, *sx0, *sx1);
    t0 = _mm_xor_si128(t0, b0);
    t1 = _mm_xor_si128(t1, b1);
    low_sqr(rl0, rl1, t0, t1);
    low_mul(&t0, &t1, *rz0, *rz1, _mm_xor_si128(ql0, ONE), ql1);
    *rl0 = _mm_xor_si128(*rl0, t0);
    *rl1 = _mm_xor_si128(*rl1, t1);

    /* Lpmq = Lpq + (Z2*X2)^2 + Zpq + (B * Z2^2)*lP */
    low_sqr(&t0, &t1, *sz0, *sz1);
    low_mul(&a0, &a1, b0, b1, t0, t1);
    low_mul(&a0, &a1, a0, a1, ql0, ql1);
    low_mul(sl0, sl1, *sx0, *sx1, *sz0, *sz1);
    low_sqr(sl0, sl1, *sl0, *sl1);
    *sl0 = _mm_xor_si128(*sl0, *rl0); *sl1 = _mm_xor_si128(*sl1, *rl1);
    *sl0 = _mm_xor_si128(*sl0, *rz0); *sl1 = _mm_xor_si128(*sl1, *rz1);
    *sl0 = _mm_xor_si128(*sl0, a0); *sl1 = _mm_xor_si128(*sl1, a1);
    /* Xpmq = Xpq + C * Z2^2 */
    low_mul(&c0, &c1, c0, c1, t0, t1);
    *sx0 = _mm_xor_si128(*rx0, c0); *sx1 = _mm_xor_si128(*rx1, c1);
    /* Zpmq = Zpq + B * Z2^2 */
    low_mul(&b0, &b1, b0, b1, t0, t1);
    *sz0 = _mm_xor_si128(*rz0, b0); *sz1 = _mm_xor_si128(*rz1, b1);

    /* end */
    return;
}

/* [pm > p] doubling-and-addition-addition */
void eca_dbl_add_add(__m128i *rx0, __m128i *rx1,
                 __m128i *rl0, __m128i *rl1,
                 __m128i *rz0, __m128i *rz1,
                 __m128i px0, __m128i px1,
                 __m128i pl0, __m128i pl1,
                 __m128i pz0, __m128i pz1,
                 __m128i qx0, __m128i qx1,
                 __m128i ql0, __m128i ql1,
                 __m128i sx0, __m128i sx1,
                 __m128i sl0, __m128i sl1) {
    /* var */
    __m128i d0, d1, e0, e1, f0, f1, s0, s1;
    __m128i t0, t1, a0, a1, b0, b1, c0, c1;
    __m128i x0, x1;
    __m128i ONE;

    /* init */
    ONE = _mm_set_epi64x(0x0, 0x1);

    /* point doubling-and-addition */
    /* D */
    low_sqr(&d0, &d1, pl0, pl1);
    /* E */
    low_sqr(&e0, &e1, pz0, pz1);
    /* F */
    f0 = _mm_xor_si128(ql0, ONE); f1 = ql1;
    /* S */
    low_mul(&s0, &s1, qx0, qx1, e0, e1);

    /* T */
    low_mul(&t0, &t1, pl0, pl1, pz0, pz1);
    t0 = _mm_xor_si128(t0, d0); t1 = _mm_xor_si128(t1, d1);
    low_mul_00u(&x0, &x1, e0, e1);
    t0 = _mm_xor_si128(t0, x0); t1 = _mm_xor_si128(t1, x1);
    /* A */
    low_sqr(&a0, &a1, px0, px1);
    low_mul(&a0, &a1, a0, a1, e0, e1);
    x0 = _mm_xor_si128(f0, _mm_slli_si128(ONE, 8)); x1 = f1;
    low_mul(&x0, &x1, x0, x1, e0, e1);
    x0 = _mm_xor_si128(x0, d0); x1 = _mm_xor_si128(x1, d1);
    low_mul(&x0, &x1, x0, x1, t0, t1);
    a0 = _mm_xor_si128(a0, x0); a1 = _mm_xor_si128(a1, x1);
    /* B */
    b0 = _mm_xor_si128(s0, t0); b1 = _mm_xor_si128(s1, t1);
    low_sqr(&b0, &b1, b0, b1);

    /* rx */
    low_sqr(rx0, rx1, a0, a1);
    low_mul(rx0, rx1, *rx0, *rx1, s0, s1);
    /* rz */
    low_mul(rz0, rz1, a0, a1, b0, b1);
    low_mul(rz0, rz1, *rz0, *rz1, e0, e1);
    /* rl */
    *rl0 = _mm_xor_si128(a0, b0); *rl1 = _mm_xor_si128(a1, b1);
    low_sqr(rl0, rl1, *rl0, *rl1);
    low_mul(rl0, rl1, *rl0, *rl1, t0, t1);
    f0 = _mm_xor_si128(sl0, f0); f1 = _mm_xor_si128(sl1, f1);
    low_mul(&x0, &x1, f0, f1, *rz0, *rz1);
    a0 = _mm_xor_si128(*rl0, x0); a1 = _mm_xor_si128(*rl1, x1);

    /* point addition */
    /* T */
    low_mul(&t0, &t1, sx0, sx1, *rz0, *rz1);

    /* A */
    /* B */
    b0 = _mm_xor_si128(*rx0, t0); b1 = _mm_xor_si128(*rx1, t1);
    low_sqr(&b0, &b1, b0, b1);
    /* C */
    low_mul(&x0, &x1, a0, a1, t0, t1);
    /* D */
    low_mul(&d0, &d1, a0, a1, b0, b1);

    /* rx */
    low_mul(&c0, &c1, a0, a1, *rx0, *rx1);
    low_mul(rx0, rx1, x0, x1, c0, c1);
    /* rz */
    low_mul(rz0, rz1, d0, d1, *rz0, *rz1);
    /* rl */
    c0 = _mm_xor_si128(c0, b0); c1 = _mm_xor_si128(c1, b1);
    low_sqr(&c0, &c1, c0, c1);
    x0 = _mm_xor_si128(sl0, ONE); x1 = sl1;
    low_mul(&x0, &x1, x0, x1, *rz0, *rz1);
    *rl0 = _mm_xor_si128(c0, x0); *rl1 = _mm_xor_si128(c1, x1);

    /* end */
    return;
}
