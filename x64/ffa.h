/* finite field arithmetic */

/* modular reduction */
/* [F_q^2] general reduction modulo X^128 + X^64 + X */
#define low_red_128_064_001(op0,op1,op2,op3,tp0,tp1,re0,re1)\
    tp0 = _mm_xor_si128(op2, op3);\
    tp1 = _mm_xor_si128(tp0, _mm_srli_epi64(op3, 63));\
    tp1 = _mm_slli_epi64(tp1, 1);\
    re0 = _mm_xor_si128(op0, tp1);\
    re1 = _mm_xor_si128(tp0, op1);\
    re1 = _mm_xor_si128(re1, _mm_slli_epi64(op3, 1));\
    re1 = _mm_xor_si128(re1, _mm_srli_epi64(op2, 63));

/* [F_q^2] mult x Fq reduction modulo X^128 + X^64 + X */
#define low_red_128_064_001_fq1(op0,op1,op2,tp0,re0,re1)\
    tp0 = _mm_slli_epi64(op2, 1);\
    re0 = _mm_xor_si128(op0, tp0);\
    re1 = _mm_xor_si128(op2, op1);\
    re1 = _mm_xor_si128(re1, _mm_srli_epi64(op2, 63));

/* [F_q^2] squaring reduction modulo X^128 + X^64 + X */
#define low_red_128_064_001_sqr(op0,op1,op2,op3,tp0,tp1,re0,re1)\
    tp0 = _mm_xor_si128(op2, op3);\
    tp1 = _mm_slli_epi64(tp0, 1);\
    re0 = _mm_xor_si128(op0, tp1);\
    re1 = _mm_xor_si128(tp0, op1);\
    re1 = _mm_xor_si128(re1, _mm_slli_epi64(op3, 1));

/* [F_q^2] reduction of a 128-bit polynomial modulo X^127 + X^63 + 1 */
#define low_red_127_063_000(op0,op1,tp0)\
    tp0 = _mm_srli_epi64(op1, 63);\
    op0 = _mm_xor_si128(op0, tp0);\
    tp0 = _mm_slli_epi64(tp0, 63);\
    op0 = _mm_xor_si128(op0, tp0);\
    op1 = _mm_xor_si128(op1, tp0);

/* [F_q  ] general reduction modulo X^128 + X^64 + X */
#define low_red_128_064_001_bas(op0,op1,tp0,tp1,tp2,re0)\
    tp0 = _mm_alignr_epi8(op1, op1, 8);\
    tp1 = _mm_xor_si128(op1, tp0);\
    tp2 = _mm_srli_epi64(op1, 63);\
    tp1 = _mm_xor_si128(tp1, tp2);\
    tp2 = _mm_unpackhi_epi64(tp1, op1);\
    tp2 = _mm_slli_epi64(tp2, 1);\
    tp1 = _mm_slli_si128(tp1, 8);\
    tp1 = _mm_xor_si128(tp1, tp2);\
    re0 = _mm_xor_si128(op0, tp1);

/* [F_q  ] squaring reduction modulo X^128 + X^64 + X */
#define low_red_128_064_001_sqr_bas(op0,op1,tp0,tp1,re0)\
    tp0 = _mm_alignr_epi8(op1, op1, 8);\
    tp0 = _mm_xor_si128(tp0, op1);\
    tp1 = _mm_unpackhi_epi64(tp0, op1);\
    tp1 = _mm_slli_epi64(tp1, 1);\
    tp0 = _mm_slli_si128(tp0, 8);\
    tp0 = _mm_xor_si128(tp0, tp1);\
    re0 = _mm_xor_si128(op0, tp0);

/* multiplication */
/* [F_q^2] karatsuba algorithm step (middle term addition is not included) */
#define low_kts_stp(op0,op1,op2,op3,op4,op5,re0,re1,re2,ord)\
    re0 = _mm_clmulepi64_si128(op0, op1, ord);\
    re1 = _mm_clmulepi64_si128(op2, op3, ord);\
    re2 = _mm_clmulepi64_si128(op4, op5, ord);\
    re1 = _mm_xor_si128(re1, re0);\
    re1 = _mm_xor_si128(re1, re2);

/* [F_q^2] Karatsuba multiplication */
void low_mul(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01,
		__m128i op10, __m128i op11) {
	/* var */
	__m128i a00, a01, a02;
	__m128i k00, k01, k02;
	__m128i r00, r01, r02, r10, r11, r12, r20, r21, r22;
	__m128i rbe, rga;
	__m128i dal, dbe, dga, dde, dt0, dt1;

	/* karatsuba: pre */
	/* high level: (b_0 + b_1) | (a_0 + a_1) (LO) > a00 */
	a00 = _mm_unpacklo_epi64(op00, op10);
	a01 = _mm_unpackhi_epi64(op00, op10);
	a00 = _mm_xor_si128(a00, a01);

	/* high level: (b_0 + b_1) | (a_0 + a_1) (HI) > a01 */
	a01 = _mm_unpacklo_epi64(op01, op11);
	a02 = _mm_unpackhi_epi64(op01, op11);
	a01 = _mm_xor_si128(a01, a02);

	/* low level: (a_0 + a_1) */
	k00 = _mm_xor_si128(op00, op01);
	k01 = _mm_xor_si128(a00, a01);
	k02 = _mm_xor_si128(op10, op11);

	/* partial karatsuba multiplication */
	low_kts_stp(op00, op10, k00, k02, op01, op11, r00, r01, r02, 0x00);	/* a0xb0 */
	low_kts_stp(a00, a00, k01, k01, a01, a01, r10, r11, r12, 0x01);	/* (a0+a1) x (b0+b1) */
	low_kts_stp(op00, op10, k00, k02, op01, op11, r20, r21, r22, 0x11);	/* a1xb1 */

	/* karatsuba: final sum (the middle term is computed separately, and then reorganized) */
	/* imaginary part */
	r10 = _mm_xor_si128(r10, r00);	/* low term */
	r11 = _mm_xor_si128(r11, r01);	/* middle term */
	r12 = _mm_xor_si128(r12, r02);	/* high term */

	/* real part */
	r00 = _mm_xor_si128(r20, r00);	/* low term */
	r01 = _mm_xor_si128(r21, r01);	/* middle term */
	r02 = _mm_xor_si128(r22, r02);	/* high term */

	rbe = _mm_unpacklo_epi64(r01, r11);
	rga = _mm_unpackhi_epi64(r01, r11);

	/* reduction: pre */
	dal = _mm_unpacklo_epi64(r00, r10);
	dbe = _mm_unpackhi_epi64(r00, r10);
	dga = _mm_unpacklo_epi64(r02, r12);
	dde = _mm_unpackhi_epi64(r02, r12);

	/* karatsuba: final sum (middle term is added to the result values) */
	dbe = _mm_xor_si128(dbe, rbe);
	dga = _mm_xor_si128(dga, rga);

	/* reduction */
	low_red_128_064_001(dal, dbe, dga, dde, dt0, dt1, *re_0, *re_1);

	/* end */
	return;
}

/* [F_q^2] multiplication by (1 + u) */
void low_mul_01u(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	/* var */
	__m128i t00, t01;

	/* multiplication */
	/* (a0 + a1u) * (1 + u) = (a0 + a1) + a0u */
	t00 = _mm_slli_si128(op00, 8);
	t01 = _mm_slli_si128(op01, 8);

	*re_0 = _mm_xor_si128(op00, t00);
	*re_1 = _mm_xor_si128(op01, t01);

	*re_0 = _mm_alignr_epi8(*re_0, *re_0, 8);
	*re_1 = _mm_alignr_epi8(*re_1, *re_1, 8);

	/* end */
	return;
}

/* [F_q^2] multiplication by (0 + u) */
void low_mul_00u(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	/* var */
	__m128i t00, t01;

	/* multiplication */
	/* (a0 + a1u) * (0 + u) = a1 + (a0 + a1)u */
	t00 = _mm_srli_si128(op00, 8);
	t01 = _mm_srli_si128(op01, 8);

	*re_0 = _mm_xor_si128(op00, t00);
	*re_1 = _mm_xor_si128(op01, t01);

	*re_0 = _mm_alignr_epi8(*re_0, *re_0, 8);
	*re_1 = _mm_alignr_epi8(*re_1, *re_1, 8);

	/* end */
	return;
}

/* [F_q^2] multiplication by (x^27 + u) */
void low_mul_027(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	/* var */
	__m128i l01, r00;
	__m128i alp, bet, gam;
	__m128i t00, t01;

	/* multiplication */
	alp = _mm_slli_epi64(op00, 27);
	l01 = _mm_slli_epi64(op01, 27);

	r00 = _mm_srli_epi64(op00, 37);
	gam = _mm_srli_epi64(op01, 37);

	bet = _mm_xor_si128(l01, r00);

	/* reduction */
	*re_1 = _mm_xor_si128(bet, gam);
	gam = _mm_slli_epi64(gam, 1);
	*re_0 = _mm_xor_si128(alp, gam);

	/* end */
	return;
}

/* [F_q^2] multiplication by (x^27 + u) */
void low_mul_27u(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	/* var */
	__m128i l01, r00;
	__m128i alp, bet, gam;
	__m128i t00, t01;

	/* multiplication */
	alp = _mm_slli_epi64(op00, 27);
	l01 = _mm_slli_epi64(op01, 27);

	r00 = _mm_srli_epi64(op00, 37);
	gam = _mm_srli_epi64(op01, 37);

	bet = _mm_xor_si128(l01, r00);

	/* reduction */
	bet = _mm_xor_si128(bet, gam);
	gam = _mm_slli_epi64(gam, 1);
	alp = _mm_xor_si128(alp, gam);

	/* final multiplication */
	t00 = _mm_srli_si128(op00, 8);
	t01 = _mm_srli_si128(op01, 8);

	t00 = _mm_xor_si128(t00, op00);
	t01 = _mm_xor_si128(t01, op01);

	t00 = _mm_alignr_epi8(t00, t00, 8);
	t01 = _mm_alignr_epi8(t01, t01, 8);

	*re_0 = _mm_xor_si128(t00, alp);
	*re_1 = _mm_xor_si128(t01, bet);

	/* end */
	return;
}

/* [F_q^2] multiplication by (b \in Fq) */
void low_mul_fq1(__m128i *c_00, __m128i *c_01, __m128i a_00, __m128i a_01,
		__m128i b_00) {
	/* var */
	__m128i re00, re01, im00, im01;
	__m128i real, rebe, rega, rede;
	__m128i tmp0;

	/* multiplication */
	re00 = _mm_clmulepi64_si128(a_00, b_00, 0x00);
	re01 = _mm_clmulepi64_si128(a_01, b_00, 0x00);
	im00 = _mm_clmulepi64_si128(a_00, b_00, 0x01);
	im01 = _mm_clmulepi64_si128(a_01, b_00, 0x01);

	/* reduction: pre */
	real = _mm_unpacklo_epi64(re00, im00);
	rebe = _mm_unpackhi_epi64(re00, im00);
	rega = _mm_unpacklo_epi64(re01, im01);
	rede = _mm_unpackhi_epi64(re01, im01);

	/* multiplication: post */
	rebe = _mm_xor_si128(rebe, rega);
	rega = rede;

	/* reduction */
	low_red_128_064_001_fq1(real, rebe, rega, tmp0, *c_00, *c_01);

	/* end */
	return;
}

/* [F_q  ] Karatsuba multiplication */
void low_mul_bas(__m128i *re_0, __m128i op00, __m128i op10) {
	/* var */
	__m128i a00, a01, a02;
	__m128i k00, k01, k02;
	__m128i sal, sbe;

	/* karatsuba: pre */
	a00 = _mm_unpacklo_epi64(op00, op10);
	a01 = _mm_unpackhi_epi64(op00, op10);
	a00 = _mm_xor_si128(a00, a01);

	/* karatsuba */
	k00 = _mm_clmulepi64_si128(op00, op10, 0x00);
	k01 = _mm_clmulepi64_si128(a00, a00, 0x10);
	k02 = _mm_clmulepi64_si128(op00, op10, 0x11);
	k01 = _mm_xor_si128(k01, k00);
	k01 = _mm_xor_si128(k01, k02);

	/* katatsuba: post */
	sal = _mm_xor_si128(k00, _mm_slli_si128(k01, 8));
	sbe = _mm_xor_si128(k02, _mm_srli_si128(k01, 8));

	/* reduction */
	low_red_128_064_001_bas(sal, sbe, a00, a01, a02, *re_0);

	/* end */
	return;
}

/* squaring */
/* [F_q^2] squaring */
void low_sqr(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	/* var */
	__m128i a00, a01;
	__m128i sal, sbe, sga, sde, st0, st1;
	__m128i dt0, dt1;

	/* pre */
	a00 = _mm_shuffle_epi32(op00, 0xD8);
	a01 = _mm_shuffle_epi32(op01, 0xD8);

	/* squaring */
	sal = _mm_clmulepi64_si128(a00, a00, 0x00);
	sbe = _mm_clmulepi64_si128(a00, a00, 0x11);
	sga = _mm_clmulepi64_si128(a01, a01, 0x00);
	sde = _mm_clmulepi64_si128(a01, a01, 0x11);

	/* reduce */
	low_red_128_064_001_sqr(sal, sbe, sga, sde, st0, st1, *re_0, *re_1);

	/* squaring: final sum */
	st0 = _mm_srli_si128(*re_0, 8);
	st1 = _mm_srli_si128(*re_1, 8);
	*re_0 = _mm_xor_si128(*re_0, st0);
	*re_1 = _mm_xor_si128(*re_1, st1);

	/* end */
	return;
}

/* [F_q  ] squaring */
void low_sqr_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i sal, sbe;
	__m128i t00, t01;

	/* squaring */
	sal = _mm_clmulepi64_si128(op00, op00, 0x00);
	sbe = _mm_clmulepi64_si128(op00, op00, 0x11);

	/* reduce */
	low_red_128_064_001_sqr_bas(sal, sbe, t00, t01, *re_0);

	/* end */
	return;
}

#define low_sqr_003_stp(fir)\
    sal = _mm_clmulepi64_si128(fir, fir, 0x00);\
    sbe = _mm_clmulepi64_si128(fir, fir, 0x11);\
    low_red_128_064_001_sqr_bas(sal,sbe,t00,t01,*re_0);\
    sal = _mm_clmulepi64_si128(*re_0, *re_0, 0x00);\
    sbe = _mm_clmulepi64_si128(*re_0, *re_0, 0x11);\
    low_red_128_064_001_sqr_bas(sal,sbe,t00,t01,*re_0);\
    sal = _mm_clmulepi64_si128(*re_0, *re_0, 0x00);\
    sbe = _mm_clmulepi64_si128(*re_0, *re_0, 0x11);\
    low_red_128_064_001_sqr_bas(sal,sbe,t00,t01,*re_0);

/* [F_q  ] three squarings */
void low_sqr_003_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i sal, sbe;
	__m128i t00, t01;

	/* squarings */
	low_sqr_003_stp(op00);

	/* end */
	return;
}

/* [F_q  ] six squarings */
void low_sqr_006_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i sal, sbe;
	__m128i t00, t01;

	/* squarings */
	low_sqr_003_stp(op00);
	low_sqr_003_stp(*re_0);		/* 006 */

	/* end */
	return;
}

/* [F_q  ] fifteen squarings */
void low_sqr_015_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i sal, sbe;
	__m128i t00, t01;

	/* squarings */
	low_sqr_003_stp(op00);
	low_sqr_003_stp(*re_0);		/* 006 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 012 */
	low_sqr_003_stp(*re_0);		/* 015 */

	/* end */
	return;
}

/* [F_q  ] thirty squarings */
void low_sqr_030_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i sal, sbe;
	__m128i t00, t01;

	/* squarings */
	low_sqr_003_stp(op00);
	low_sqr_003_stp(*re_0);		/* 006 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 012 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 018 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 024 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 030 */

	/* end */
	return;
}

/* [F_q  ] sixty-three squarings */
void low_sqr_063_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i sal, sbe;
	__m128i t00, t01;

	/* squarings */
	low_sqr_003_stp(op00);
	low_sqr_003_stp(*re_0);		/* 006 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 012 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 018 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 024 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 030 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 036 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 042 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 048 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 054 */
	low_sqr_003_stp(*re_0);
	low_sqr_003_stp(*re_0);		/* 060 */
	low_sqr_003_stp(*re_0);		/* 063 */

	/* end */
	return;
}

#include "sqr_tbl.inc"

//MULTI SQUARING F_2^127 (6x)
void low_sqr06(__m128i *b, __m128i _a) {
	__m128i r0;
	uint64_t *p, a[2];
	int i;

    _mm_store_si128((__m128i *)a, _a);
	r0 = _mm_setzero_si128();

	for (i=0;i<16;i++) {
		p = tbl_sqr06[i][(a[0]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));

		p = tbl_sqr06[i+16][(a[1]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
	}

    *b = r0;
}

//MULTI SQUARING F_2^127 (12x)
void low_sqr12(__m128i *b, __m128i _a) {
    __m128i r0;
	uint64_t *p, a[2];
	int i;

    _mm_store_si128((__m128i *)a, _a);
	r0 = _mm_setzero_si128();

	for (i=0;i<16;i++) {
		p = tbl_sqr12[i][(a[0]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));

		p = tbl_sqr12[i+16][(a[1]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
	}

    *b = r0;
}

//MULTI SQUARING F_2^127 (24x)
void low_sqr24(__m128i *b, __m128i _a) {
    __m128i r0;
	uint64_t *p, a[2];
	int i;

    _mm_store_si128((__m128i *)a, _a);
	r0 = _mm_setzero_si128();

	for (i=0;i<16;i++) {
		p = tbl_sqr24[i][(a[0]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));

		p = tbl_sqr24[i+16][(a[1]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
	}

    *b = r0;
}

//MULTI SQUARING F_2^127 (48x)
void low_sqr48(__m128i *b, __m128i _a) {
    __m128i r0;
	uint64_t *p, a[2];
	int i;

    _mm_store_si128((__m128i *)a, _a);
	r0 = _mm_setzero_si128();

	for (i=0;i<16;i++) {
		p = tbl_sqr48[i][(a[0]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));

		p = tbl_sqr48[i+16][(a[1]>>(4*i))&0x0F];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
	}

	 *b = r0;
}

void low_inv_tbl(__m128i *_b, __m128i a) {
	__m128i tmp, b, a2_6, a2_24;

	//a^(2^2-2) * a = a^(2^2-1)
	low_sqr_bas(&b, a);
	low_mul_bas(&b, b, a);

	//a^(2^3-2) * a = a^(2^3-1)
	low_sqr_bas(&b, b);
	low_mul_bas(&b, b, a);

	//a^(2^6-2^3) * a^(2^3-1) = a^(2^6-1)
    low_sqr_003_bas(&a2_6, b);	/* 003 > 006 */
	low_mul_bas(&a2_6, a2_6, b);

	//a^(2^12-2^6) * a^(2^6-1) = a^(2^12-1)
	//low_sqi(b, a2_6, 6);
	low_sqr06(&b, a2_6);
	low_mul_bas(&b, b, a2_6);

	//a^(2^24-2^12) * a^(2^12-1) = a^(2^24-1)
	low_sqr12(&a2_24, b);
	low_mul_bas(&a2_24, a2_24, b);

	//a^(2^48-2^24) * a^(2^24-1) = a^(2^48-1)
	low_sqr24(&b, a2_24);
	low_mul_bas(&b, b, a2_24);

	//a^(2^96-2^48) * a^(2^48-1) = a^(2^96-1)
	low_sqr48(&tmp, b);
	low_mul_bas(&b, tmp, b);

	//a^(2^120-2^24) * a^(2^24-1) = a^(2^120-1)
	low_sqr24(&b, b);
	low_mul_bas(&b, b, a2_24);

	//a^(2^126-2^6) * a^(2^6-1) = a^(2^126-1)
	//low_sqi(b, b, 6);
	low_sqr06(&b, b);
	low_mul_bas(&b, b, a2_6);

	//a^(2^127-2)
	low_sqr_bas(_b, b);
}

/* [F_q^2] inversion */
void low_inv_var(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	/* var */
	__m128i are, aim, cre, cim;
	__m128i t00, t01, t02;

	/* inversion: pre */
	are = _mm_unpacklo_epi64(op00, op01);
	aim = _mm_unpackhi_epi64(op00, op01);

	/* inversion */
	t00 = _mm_xor_si128(are, aim);	/* t00 = a_r + a_i */
	low_mul_bas(&t01, are, aim);	/* t01 = a_r * a_i */
	low_sqr_bas(&t02, t00);		/* t02 = (a_r + a_i)^2 */
	t01 = _mm_xor_si128(t01, t02);	/* t = a_r * a_i + (a_r + a_i)^2 */

	low_inv_tbl(&t01, t01);		/* t = t^-1 */

	low_mul_bas(&cre, t00, t01);	/* c_r = (a_r + a_i) * t */
	low_mul_bas(&cim, aim, t01);	/* c_i = a_i * t */

	/* inversion: post */
	*re_0 = _mm_unpacklo_epi64(cre, cim);
	*re_1 = _mm_unpackhi_epi64(cre, cim);

	/* end */
	return;
}

/* inversion */
/* [F_q  ] inversion */
void low_inv_bas(__m128i *re_0, __m128i op00) {
	/* var */
	__m128i a00, a01, a03;

	/* itoh-tsujii: 1-2-3-6-12-15-30-60-63-126 */
	low_sqr_bas(&a01, op00);	/* 001 > 002 */
	low_mul_bas(&a00, a01, op00);

	low_sqr_bas(&a01, a00);		/* 002 > 003 */
	low_mul_bas(&a03, a01, op00);

	low_sqr_003_bas(&a01, a03);	/* 003 > 006 */
	low_mul_bas(&a00, a01, a03);

	low_sqr_006_bas(&a01, a00);	/* 006 > 012 */
	low_mul_bas(&a00, a01, a00);

	low_sqr_003_bas(&a01, a00);	/* 012 > 015 */
	low_mul_bas(&a00, a01, a03);

	low_sqr_015_bas(&a01, a00);	/* 015 > 030 */
	low_mul_bas(&a00, a01, a00);

	low_sqr_030_bas(&a01, a00);	/* 030 > 060 */
	low_mul_bas(&a00, a01, a00);

	low_sqr_003_bas(&a01, a00);	/* 060 > 063 */
	low_mul_bas(&a00, a01, a03);

	low_sqr_063_bas(&a01, a00);	/* 063 > 126 */
	low_mul_bas(&a00, a01, a00);

	low_sqr_bas(re_0, a00);		/* 126 > 127 */

	/* end */
	return;
}

void _low_inv_gcd(__m128i *re_0, __m128i op00) {
    uint64_t c0, c1, f[2], v[2];
	int64_t delta = -1;
    uint64_t g0 = 0x8000000000000001, g1 = 0x8000000000000000;
	uint64_t u0 = 1, u1 = 0, v0 = 0, v1 = 0;

	_mm_storeu_si128((__m128i *)f, op00);
	uint64_t f0 = f[0];
	uint64_t f1 = f[1];

	for (int i = 0; i <= 2*127; i++) {
        c0 = delta >> 63;
        c1 = -(f0 & 1);

		f0 ^= (g0 & c1);
		f1 ^= (g1 & c1);

		u0 ^= (v0 & c1);
		u1 ^= (v1 & c1);

		g0 ^= (f0 & c0);
		g1 ^= (f1 & c0);

		v0 ^= (u0 & c0);
		v1 ^= (u1 & c0);

        c0 &= c1;
        delta = (delta ^ c0) - 1;

		f0 = (f0 >> 1) | (f1 << 63);
		f1 = (f1 >> 1);

		uint64_t _u = u0 & 1;
		u0 ^= _u;
		u0 ^= _u << 63;
		u1 ^= _u << 63;

		u0 = (u0 >> 1) | (u1 << 63);
		u1 = (u1 >> 1);
	}

    v[0] = v0;
    v[1] = v1;
	*re_0 = _mm_loadu_si128((__m128i *)v);
	//gf_out(*re_0);
}

static inline ulong revbin(ulong x) {
    int i;

    union rev_t {
        uint8_t bytes[8];
        uint64_t value;
    } rev;

    rev.value = x;
    for (i = 0; i < 8; i++)
        rev.bytes[i] = (rev.bytes[i] * 0x0202020202ULL & 0x010884422010ULL) % 1023;
    return __builtin_bswap64(rev.value);
}

void low_inv_gcd(__m128i *re_0, __m128i op00) {
    uint64_t c[2] = {0x2000000000000000, 0x3000000000000000};
    uint64_t f[2] = {0x8000000000000001, 0x8000000000000000};
    uint64_t c0, c1, g[2], v[2];
	uint64_t delta = 1, mask = (uint64_t)1 << 63;
    uint64_t t1, t0, u0 = 1, u1 = 0, v0 = 0, v1 = 0;
    uint64_t f0 = 0x1;
	uint64_t f1 = 0x8000000000000001;

    /* Store and fully reduce in the base field. */
	_mm_storeu_si128((__m128i *)g, op00);
    g[0] ^= g[1] & mask;
    g[0] ^= g[1] >> 63;
    g[1] ^= g[1] & mask;

    uint64_t g0 = revbin(g[1]);
	uint64_t g1 = revbin(g[0]);
    g0 = (g0 >> 1) | (g1 << 63);
    g1 = (g1 >> 1);

	for (int i = 2*127 - 1; i > 0; i--) {
        /* If delta > 0 && g0 & 1. */
        c0 = (delta >> 63);
        c1 = (g0 & 1);
        c0 = (c0 ^ 1) & c1;

        delta = (delta ^ -c0) + c0 + 1;

        t0 = (f0 ^ g0) & -c0;
        f0 ^= t0;
        g0 ^= t0;
        t1 = (f1 ^ g1) & -c0;
        f1 ^= t1;
        g1 ^= t1;

        t0 = (u0 ^ v0) & -c0;
        u0 ^= t0;
        v0 ^= t0;
        t1 = (u1 ^ v1) & -c0;
        u1 ^= t1;
        v1 ^= t1;

        c0 = -(f0 & 1);
        c1 = -(g0 & 1);

        g0 = (g0 & c0) ^ (f0 & c1);
        g1 = (g1 & c0) ^ (f1 & c1);
        u0 = (u0 & c0) ^ (v0 & c1);
        u1 = (u1 & c0) ^ (v1 & c1);

        g0 = (g0 >> 1) | (g1 << 63);
		g1 = (g1 >> 1);

        v1 = (v1 << 1) | (v0 >> 63);
        v0 = (v0 << 1);
        v0 ^= v1 & mask;
        v0 ^= v1 >> 63;
        v1 ^= v1 & mask;
	}

    c0 = v0 & 1;
    // Divite by x. */
    v[0] = (v0 >> 1) | (v1 << 63);
    v[1] = (v1 >> 1);
    // Revert order and handle carry out of the shift above.
    v1 = revbin(v[0]);
	v0 = revbin(v[1]);
    v[0] = (v0 >> 1) | (v1 << 63);
    v[1] = (v1 >> 1);
    v[0] ^= c0;
    v[1] ^= c0;

	*re_0 = _mm_loadu_si128((__m128i *)v);
}

void _low_inv_vec(__m128i *re_0, __m128i op00) {
    uint64_t f0, d0, f[2], u[2] = {1,0}, v[2] = {0}, s[2] = {0x8000000000000001,0x8000000000000000};
    int32_t delta = -1;

    _mm_storeu_si128((__m128i *)f, op00);
    __m128i _s = _mm_set_epi64x(s[1], s[0]);
    __m128i _v = _mm_set_epi64x(v[1], v[0]);
    __m128i _1 = _mm_set_epi64x(0x0, 0x1);
    __m128i _2 = _mm_set_epi64x(0x1, 0x1);
    __m128i _f = op00;
    __m128i _u = _1;
    __m128i _m = _mm_set_epi32(delta, delta, delta, delta);
    __m128i _d0 = _m;

    for (int i = 0; i <= 2*127; i++) {
        __m128i _f0, t;

        //_f0 = _mm_set_epi64x(-f0, -f0);
        //STR((__m128i *)(f+0), _f);
        t = _mm_and_si128(_f, _1);
        t = _mm_xor_si128(t, _mm_slli_si128(t, 8));
        _f0 = _mm_cmpeq_epi64(t, _2);
        //f0 = f[0] & 1; d0 = delta >> 31;

        //_d0 = _mm_set_epi32(delta, delta, delta, delta);
        t = _mm_srai_epi32(_d0, 31);
        //printf("%X\n", d0);
        //uint32_t _t[4];
        //STR((__m128i *)(_t+0), _d0);
        //printf("%X %X\n", t[0], t[1]);
        //exit(0);

        //f[0] ^= (s[0] & -f0);
        //f[1] ^= (s[1] & -f0);
        _f = _mm_xor_si128(_f, _mm_and_si128(_s, _f0));

        //u[0] ^= (v[0] & -f0);
        //u[1] ^= (v[1] & -f0);
        _u = _mm_xor_si128(_u, _mm_and_si128(_v, _f0));

        //s[0] ^= (f[0] & d0);
        //s[1] ^= (f[1] & d0);
        _s = _mm_xor_si128(_s, _mm_and_si128(_f, t));

        //v[0] ^= (u[0] & d0);
        //v[1] ^= (u[1] & d0);
        _v = _mm_xor_si128(_v, _mm_and_si128(_u, t));

        //printf("%d\n", delta);
        //STR((__m128i *)(_t+0), _d0);
        //printf("%d %d %d %d\n", _t[3], _t[2], _t[1], _t[0]);
        //delta = SELECT(delta, -delta, f0 & -d0);
        //_d0 = _mm_sign_epi32(_d0, _mm_and_si128(_f0, t));
        _d0 = _mm_blendv_epi8(_d0, _mm_sign_epi32(_d0, _m), _mm_and_si128(_f0, t));
        //printf("%d\n", delta);
        //STR((__m128i *)(_t+0), _d0);
        //printf("%d %d %d %d\n", _t[3], _t[2], _t[1], _t[0]);

        //STR((__m128i *)(u+0), _u);

        //f[0] = (f[0] >> 1) | (f[1] << 63);
        //f[1] = (f[1] >> 1);
        _f = _mm_xor_si128(_mm_srli_epi64(_f, 1), _mm_slli_epi64(_mm_srli_si128(_f, 8), 63));

        //uint64_t u0 = u[0] & 1;
        //u[0] ^= u0;
        //u[0] ^= u0 << 63;
        //u[1] ^= u0 << 63;
        t = _mm_and_si128(_u, _1);
        _u = _mm_xor_si128(_u, t);
        t = _mm_slli_epi64(t, 63);
        _u = _mm_xor_si128(_u, t);
        t = _mm_slli_si128(t, 8);
        _u = _mm_xor_si128(_u, t);

        //u[0] = (u[0] >> 1) | (u[1] << 63);
        //u[1] = (u[1] >> 1);
        _u = _mm_xor_si128(_mm_srli_epi64(_u, 1), _mm_slli_epi64(_mm_srli_si128(_u, 8), 63));

        //delta--;
        _d0 = _mm_add_epi32(_d0, _m);
    }

    *re_0 = _v;//_mm_loadu_si128((__m128i *)v);
    //gf_out(*re_0);
}

#define CTIME_IF(mask,then,else)  ((mask&(then)) | (~mask&(else) ))

#define CSWAP(U, V, T, C)      (T = ((U) ^ (V)) & -C, U ^= T, V^= T)

static int divstepsx(int n, uint64_t delta, uint64_t f64, uint64_t g64, uint64_t * p00, uint64_t * p01, uint64_t * p10, uint64_t * p11) {
    uint64_t u, v, q, r;
    uint64_t g0, f0;

    u = ((uint64_t) 1) << n;
    v = 0;
    q = 0;
    r = ((uint64_t) 1) << n;
    uint64_t tmp,tmp2;

    while (n > 0) {
        int64_t swap_mask = ((delta > 0) & ((g64 & 1) != 0));
        swap_mask = (swap_mask << (8-1)) >> (8-1);
        delta = (delta ^ -swap_mask) + swap_mask;
        //delta = CTIME_IF(swap_mask,-delta,delta);
        CSWAP(f64, g64, tmp, swap_mask);
        CSWAP(q, u, tmp, swap_mask);
        CSWAP(r, v, tmp, swap_mask);

        delta++;
        g0 = (uint64_t)0 - (g64 & (uint64_t) 0x1);
        f0 = (uint64_t)0 - (f64 & (uint64_t) 0x1);

        q =   (f0 & q) ^ (g0 & u);
        r =   (f0 & r) ^ (g0 & v);
        g64 = (f0 & g64) ^ (g0 & f64);
        g64 >>= 1;
        q >>= 1;
        r >>= 1;
        n--;
    } //end while
    *p00 = u;
    *p01 = v;
    *p10 = q;
    *p11 = r;

    return delta;
}


static inline __m128i right_shift_128(__m128i in) {

    __m128i a,b;
    a = _mm_srli_epi64(in,1);
    b = _mm_slli_epi64(in,8-1);

    //set the high part of b = 0
    b = _mm_unpacklo_epi64( _mm_setzero_si128(),b ); //o _mm_unpacklo_epi64???
    a = _mm_or_si128(a,b);

    return a;
}

#define CTIME_IF_128(mask,then,else)  _mm_or_si128(_mm_and_si128(mask, then) ,_mm_andnot_si128(mask, else))

static inline int divstepsx_128(int n, int delta, uint64_t f[], uint64_t g[], uint64_t *p00, uint64_t *p01, uint64_t *p10, uint64_t *p11) {

    if(n<64)
        return delta = divstepsx (n, delta, f[0], g[0], p00, p01, p10, p11);

    __m128i g0, f0, g128, f128;
    __m128i one_128 = _mm_set_epi64x((uint64_t) 1, (uint64_t) 0);
    __m128i mask_128 = _mm_set_epi64x((uint64_t) 1, (uint64_t) 1);
    __m128i zero_128 = _mm_setzero_si128();

    __m128i u, v, q, r;

    g128 = _mm_lddqu_si128((__m128i *)g);
    f128 = _mm_lddqu_si128((__m128i *)f);

    uint64_t  temp = ((uint64_t) 1)<< (n-64);
    u = _mm_set_epi64x( (uint64_t) 0, temp);
    r = _mm_set_epi64x( (uint64_t) 0, temp);
    v = _mm_setzero_si128();
    q = _mm_setzero_si128();

    __m128i delta_128 = _mm_set_epi64x((int64_t) delta,(int64_t) delta);

    __m128i tmp,tmp2;

    while (n > 0) {
        __m128i delta_mask = _mm_cmpgt_epi64(delta_128, zero_128);

        //something like [xxx....xxx | FFF....FFF] where x is the actual mask
        __m128i g128_mask = _mm_cmpeq_epi64(_mm_and_si128(g128, one_128), one_128);

        __m128i swap_mask = _mm_and_si128(delta_mask, (__m128i)_mm_shuffle_pd((__m128d) g128_mask, (__m128d) g128_mask, 3));


        delta_128 = _mm_add_epi64(_mm_xor_si128(delta_128, swap_mask), _mm_and_si128(mask_128, swap_mask));

        // delta = CTIME_IF(swap,-delta,delta);
        tmp  = CTIME_IF_128(swap_mask,g128,f128);
        tmp2 = CTIME_IF_128(swap_mask,f128,g128);


        f128  = tmp;
        g128  = tmp2;

        tmp  = CTIME_IF_128(swap_mask,q,u);
        tmp2 = CTIME_IF_128(swap_mask,u,q);


        u  = tmp;
        q  = tmp2;

        tmp  = CTIME_IF_128(swap_mask,r,v);
        tmp2 = CTIME_IF_128(swap_mask,v,r);

        v  = tmp;
        r  = tmp2;

        //delta++;
        delta_128 = _mm_add_epi64(delta_128, mask_128);

        g0 = _mm_cmpeq_epi64(_mm_and_si128(g128, one_128), one_128);
        g0 = (__m128i)_mm_shuffle_pd((__m128d) g0, (__m128d) g0, 3);

        f0 = _mm_cmpeq_epi64(_mm_and_si128(f128, one_128), one_128);
        f0 = (__m128i)_mm_shuffle_pd((__m128d) f0, (__m128d) f0, 3);

        q =   _mm_xor_si128(_mm_and_si128(f0,q), _mm_and_si128(g0,u)); //(f0 & q) ^ (g0 & u);
        r =   _mm_xor_si128(_mm_and_si128(f0,r), _mm_and_si128(g0,v)); //(f0 & r) ^ (g0 & v);
        g128 = _mm_xor_si128(_mm_and_si128(f0,g128), _mm_and_si128(g0,f128)); //(f0 & g64) ^ (g0 & f64);

        g128 = right_shift_128(g128);
        q = right_shift_128(q);
        r = right_shift_128(r);
        n--;
    } //end while

    _mm_storeu_si128((__m128i *) p00, u);
    _mm_storeu_si128((__m128i *) p01, v);
    _mm_storeu_si128((__m128i *) p10, q);
    _mm_storeu_si128((__m128i *) p11, r);

    return 0;//_mm_cvtsi128_si64x(delta_128);//_mm_extract_epi64(delta_128,1);
}

static inline __m256i right_shift_256(__m256i in) {
    __m256i a,b,c;
    a = _mm256_srli_epi64(in,1);

    b = _mm256_slli_epi64(in,8-1);

    c = _mm256_permute4x64_epi64(b,0x93);
    c = _mm256_insert_epi64(c, (uint64_t) 0, 0);

    a = _mm256_or_si256(a,c);

    return a;
}

#define CTIME_IF_256(mask,then,else)  _mm256_blendv_epi8(else, then, mask)

static inline int divstepsx_256(int n, int delta, uint64_t f[], uint64_t g[], uint64_t *p00, uint64_t *p01, uint64_t *p10, uint64_t *p11) {

    if(n<128)
        return delta = divstepsx_128 (n,delta, f, g, p00, p01, p10, p11);

    __m256i g0, f0, g256, f256;
    __m256i one_256 = _mm256_set_epi64x((uint64_t) 1, (uint64_t) 0,(uint64_t) 0, (uint64_t) 0);
    __m256i mask_256 = _mm256_set_epi64x((uint64_t) 1, (uint64_t) 1,(uint64_t) 1, (uint64_t) 1);
    __m256i zero_256 = _mm256_setzero_si256();


    __m256i u, v, q, r;

    g256 = _mm256_lddqu_si256((__m256i *)g);
    f256 = _mm256_lddqu_si256((__m256i *)f);

    uint64_t  temp = ((uint64_t) 1)<< (n-192);
    u = _mm256_set_epi64x( (uint64_t) 0, (uint64_t) 0, (uint64_t) 0, temp);
    r = _mm256_set_epi64x( (uint64_t) 0, (uint64_t) 0, (uint64_t) 0, temp);
    v = _mm256_setzero_si256();
    q = _mm256_setzero_si256();

    __m256i delta_256 = _mm256_set_epi64x((int64_t) delta,(int64_t) delta, (int64_t) delta,(int64_t) delta);

    __m256i tmp,tmp2;

    while (n > 0) {
        __m256i delta_mask = _mm256_cmpgt_epi64(delta_256, zero_256);
        __m256i g256_mask = _mm256_cmpeq_epi64(_mm256_and_si256(g256, one_256), one_256);

        __m256i swap_mask = _mm256_and_si256(delta_mask, _mm256_permute4x64_epi64(g256_mask, 0xFF));


        //need to add 1 when changing sign with the xor
        delta_256 = _mm256_add_epi64(_mm256_xor_si256(delta_256, swap_mask), _mm256_and_si256(mask_256, swap_mask));
        tmp  = CTIME_IF_256(swap_mask,g256,f256);
        tmp2 = CTIME_IF_256(swap_mask,f256,g256);

        f256  = tmp;
        g256  = tmp2;

        tmp  = CTIME_IF_256(swap_mask,q,u);
        tmp2 = CTIME_IF_256(swap_mask,u,q);


        u  = tmp;
        q  = tmp2;

        tmp  = CTIME_IF_256(swap_mask,r,v);
        tmp2 = CTIME_IF_256(swap_mask,v,r);

        v  = tmp;
        r  = tmp2;

        delta_256 = _mm256_add_epi64(delta_256, mask_256);

        __m256i maskgf_tmp = _mm256_cmpeq_epi64(_mm256_and_si256(g256, one_256), one_256);
        g0 = _mm256_permute4x64_epi64(maskgf_tmp, 0xFF);

        maskgf_tmp = _mm256_cmpeq_epi64(_mm256_and_si256(f256, one_256), one_256);
        f0 = _mm256_permute4x64_epi64(maskgf_tmp, 0xFF);

        q =   _mm256_xor_si256(_mm256_and_si256(f0,q), _mm256_and_si256(g0,u));
        r =   _mm256_xor_si256(_mm256_and_si256(f0,r), _mm256_and_si256(g0,v));
        g256 = _mm256_xor_si256(_mm256_and_si256(f0,g256), _mm256_and_si256(g0,f256));

        g256 = right_shift_256(g256);
        q = right_shift_256(q);
        r = right_shift_256(r);
        n--;

    } //end while

    _mm256_storeu_si256((__m256i *) p00, u);
    _mm256_storeu_si256((__m256i *) p01, v);
    _mm256_storeu_si256((__m256i *) p10, q);
    _mm256_storeu_si256((__m256i *) p11, r);

    return _mm256_extract_epi64(delta_256, 0);
}

/* [F_q^2] inversion */
void low_inv(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
    /* var */
    __m128i are, aim, cre, cim;
    __m128i t00, t01, t02;

    /* inversion: pre */
    //low_red_127_063_000(op00, op01, t00);
    are = _mm_unpacklo_epi64(op00, op01);
    aim = _mm_unpackhi_epi64(op00, op01);

    /* inversion */
    t00 = _mm_xor_si128(are, aim);  /* t00 = a_r + a_i */
    low_mul_bas(&t01, are, aim);    /* t01 = a_r * a_i */
    low_sqr_bas(&t02, t00);         /* t02 = (a_r + a_i)^2 */
    t01 = _mm_xor_si128(t01, t02);  /* t = a_r * a_i + (a_r + a_i)^2 */

    low_inv_bas(&t01, t01);         /* t = t^-1 */

    low_mul_bas(&cre, t00, t01);    /* c_r = (a_r + a_i) * t */
    low_mul_bas(&cim, aim, t01);    /* c_i = a_i * t */

    /* inversion: post */
    *re_0 = _mm_unpacklo_epi64(cre, cim);
    *re_1 = _mm_unpackhi_epi64(cre, cim);

    /* end */
    return;
}

void low_inv_sim(__m128i *o0, __m128i *o1, __m128i *i0, __m128i *i1, uint64_t len, int ct) {
    __m128i u0, u1, c0[len], c1[len];

	c0[0] = i0[0];
    c1[0] = i1[0];
	for(int i = 1; i < len; i++) {
        low_mul(&c0[i], &c1[i], c0[i-1], c1[i-1], i0[i], i1[i]);
	}
    if (ct) {
        low_inv(&u0, &u1, c0[len - 1], c1[len - 1]);
    } else {
        low_inv_var(&u0, &u1, c0[len - 1], c1[len - 1]);
    }
	for(int i = len-1; i >= 1; i--) {
        low_mul(&o0[i], &o1[i], u0, u1, c0[i-1], c1[i-1]);
		low_mul(&u0, &u1, u0, u1, i0[i], i1[i]);
	}
    o0[0] = u0;
    o1[0] = u1;
}

void low_sqrt(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	__m128i a0, aL, aH, uu0, uu1, vv0, vv1;
	__m128i ma, mb;
	uint64_t uu[2], vv[2];

	__m128i perm =
			_mm_set_epi32(0x0F0D0B09, 0x07050301, 0x0E0C0A08, 0x06040200);
	__m128i sqrtL =
			_mm_set_epi32(0x33322322, 0x31302120, 0x13120302, 0x11100100);
	__m128i sqrtH =
			_mm_set_epi32(0xCCC88C88, 0xC4C08480, 0x4C480C08, 0x44400400);
	__m128i maskL =
			_mm_set_epi32(0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F);
	__m128i maskH =
			_mm_set_epi32(0xF0F0F0F0, 0xF0F0F0F0, 0xF0F0F0F0, 0xF0F0F0F0);

	/* inversion: pre */
	ma = _mm_unpacklo_epi64(op00, op01);
	mb = _mm_unpackhi_epi64(op00, op01);
	ma = _mm_xor_si128(ma, mb);

	//Extraction of even (ae0) and odd (ao0) bits
	uu0 = _mm_shuffle_epi8(ma, perm);

	uu1 = _mm_and_si128(uu0, maskL);
	vv1 = _mm_and_si128(uu0, maskH);
	vv1 = _mm_srli_epi64(vv1, 4);

	uu1 = _mm_shuffle_epi8(sqrtL, uu1);
	vv1 = _mm_shuffle_epi8(sqrtH, vv1);

	uu0 = _mm_xor_si128(uu1, vv1);

	uu1 = _mm_and_si128(uu0, maskL);
	vv1 = _mm_and_si128(uu0, maskH);

	//Extraction of even (ae0) and odd (ao0) bits
	a0 = _mm_shuffle_epi8(mb, perm);

	aL = _mm_and_si128(a0, maskL);
	aH = _mm_and_si128(a0, maskH);
	aH = _mm_srli_epi64(aH, 4);

	aL = _mm_shuffle_epi8(sqrtL, aL);
	aH = _mm_shuffle_epi8(sqrtH, aH);

	a0 = _mm_xor_si128(aL, aH);

	aL = _mm_and_si128(a0, maskL);
	aH = _mm_and_si128(a0, maskH);

	//Multiplication of odd vector to constant value sqrt(x)
	//sqrt(x) = x^64 + x^32
	uu0 = _mm_unpacklo_epi64(uu1, aL);
	uu1 = _mm_unpackhi_epi64(uu1, aL);
	vv0 = _mm_unpacklo_epi64(vv1, aH);
	vv1 = _mm_unpackhi_epi64(vv1, aH);

	uu1 = _mm_slli_epi64(uu1, 4);
	vv0 = _mm_srli_epi64(vv0, 4);

	uu0 = _mm_xor_si128(uu0, uu1);
	vv0 = _mm_xor_si128(vv0, vv1);

	uu0 = _mm_xor_si128(uu0, _mm_slli_epi64(vv0, 32));	//b2b0
	vv0 = _mm_xor_si128(vv0, _mm_srli_epi64(vv0, 32));	//b3b1

	*re_0 = _mm_unpacklo_epi64(uu0, vv0);
	*re_1 = _mm_unpackhi_epi64(uu0, vv0);
}

#include "htr_tbl.inc"

/* half-trace in F_{q} */
void low_htr_tbl(__m128i *b, __m128i _a) {
	int i;
	uint64_t *p, a[2];
	uint8_t *tmp0, *tmp1;
	uint64_t a0_h, tmp;
    __m128i r0;

    _mm_store_si128((__m128i *)a, _a);
	a0_h = a[0] >> 32;
	tmp = _pdep_u64(a0_h, 0x5555555555555555);
	tmp = a[1] ^ tmp;

	//accumulator initialization
	a0_h = a0_h << 32;
	r0 = _mm_set_epi64x(0, a0_h);

	//window-8 pointer
	tmp0 = (uint8_t *) &a[0];
	tmp1 = (uint8_t *) &tmp;

	//look-up table
	for (i=0;i<4;i++) {
		p = tbl_htr[i][*tmp0++];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
		p = tbl_htr[i+8][*tmp1++];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
	}

	for (i=4;i<8;i++) {
		p = tbl_htr[i+8][*tmp1++];
		r0 = _mm_xor_si128(r0, *(__m128i *)(p));
	}

	*b = r0;
}

void low_htr(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
	__m128i ma, a0, a1, b0, b1;

	a0 = _mm_unpacklo_epi64(op00, op01);
	a1 = _mm_unpackhi_epi64(op00, op01);

	low_htr_tbl(&b1, a1);
	ma = _mm_xor_si128(a0, b1);
	ma = _mm_xor_si128(ma, a1);
	low_htr_tbl(&b0, ma);
	b1 = _mm_xor_si128(b1, _mm_and_si128(ma, _mm_set_epi64x(0, 1)));

	*re_0 = _mm_unpacklo_epi64(b0, b1);
	*re_1 = _mm_unpackhi_epi64(b0, b1);
}

void low_htr_bas(__m128i *b, __m128i a) {
	int i, j;
	__m128i t, _b;

	t = a;
	_b = a;

	for (i = 1; i <= (127 - 1) / 2; i++) {
		low_sqr_bas(&_b, _b);
		low_sqr_bas(&_b, _b);
		_b = _mm_xor_si128(_b, t);
	}
	*b = _b;
}

void low_htr_const(__m128i *re_0, __m128i *re_1, __m128i op00, __m128i op01) {
    __m128i ma, a0, a1, b0, b1;

	a0 = _mm_unpacklo_epi64(op00, op01);
	a1 = _mm_unpackhi_epi64(op00, op01);

	low_htr_bas(&b1, a1);
	ma = _mm_xor_si128(a0, b1);
	ma = _mm_xor_si128(ma, a1);
	low_htr_bas(&b0, ma);
	b1 = _mm_xor_si128(b1, _mm_and_si128(ma, _mm_set_epi64x(0, 1)));

	*re_0 = _mm_unpacklo_epi64(b0, b1);
	*re_1 = _mm_unpackhi_epi64(b0, b1);
}

void low_out(__m128i a0, __m128i a1) {
    uint64_t t[4];
    __m128i r;

    low_red_127_063_000(a0, a1, r);
    _mm_store_si128((__m128i *)(t+0), a0);
    _mm_store_si128((__m128i *)(t+2), a1);

    printf("0x%.16lX%.16lX, ", t[2], t[0]);
    printf("0x%.16lX%.16lX\n", t[3], t[1]);
}
