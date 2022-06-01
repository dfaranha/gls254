#define CEIL(A, B)	(((A) - 1) / (B) + 1)
#define MASK(B)     (((uint64_t)1 << (B)) - 1)

/* regular recoding of a 127-bit integer */
void smu_reg_rec(int8_t *dig, uint64_t *k, uint64_t w) {
	/* var */
	int64_t ki, sig, i, l = CEIL(128, w - 1);
    uint64_t mask = MASK(w);

	/* main loop */
	for (i = 0; i < l - 1; i++) {
		ki = (k[0] & mask) - ((int64_t)1 << (w - 1));
		dig[i] = ki;
        k[0] -= ki;

		k[0] = k[0] >> (w - 1);
		k[0] = k[0] | (k[1] << (64 - (w - 1)));
		k[1] = k[1] >> (w - 1);
	}
	dig[l - 1] = k[0] & mask;
}

#include <assert.h>

/* GLS recoding in constant time, using new approach. */
void gls_recoding(uint64_t k[], uint64_t k1[], uint64_t k2[], uint64_t *s1, uint64_t *s2) {
    //c_i <= (q+1+|t|)/((q-1)^2 + t^2)*2^256 <= q+1+2sqrt(q)/(q+1)^2 * 2^256 <= 2^130-1
	//So the c_i will at most be 130 bits
    const uint64_t c1[4] = {0xE334618602D4CB44, 0x0000000000000004, 0x2, 0x0};
    const uint64_t c2[4] = {0x1CCB9E79FD2B34AC, 0xFFFFFFFFFFFFFFFB, 0x1, 0x0};
    //|a_i| <= q+1+|t| <= q+sqrt(q)+1 <= 2q-1
	//So the a_i will at most be 128 bits
    const uint64_t a1[2] = {0x38CD186180B532D2, 0x8000000000000001};
    const uint64_t a2[2] = {0xC732E79E7F4ACD2C, 0x7FFFFFFFFFFFFFFE};
    uint64_t v1[2][2], v2[2][2];
    uint64_t b1[2], b2[2], p1, p2, tmp[8];

    //Step 2: Initializing basis
    /* v1 = (alpha1//2, alpha2//2), v2 = (alpha2//2, -alpha1//2) */
	//Need to remember that v2[1] is negative
    v2[1][0] = v1[0][0] = (a1[0] >> 1) | a1[1] << 63;
    v2[1][1] = v1[0][1] = (a1[1] >> 1);
    v2[0][0] = v1[1][0] = (a2[0] >> 1) | a2[1] << 63;
    v2[0][1] = v1[1][1] = (a2[1] >> 1);

    /* Step 3: Computing b1 = (c1*k) // 2^d, b2 = (c2*k) // 2^d for d = 256. */
    bn_muln_low(tmp, c1, k, 4);
    b1[0] = tmp[4];
    b1[1] = tmp[5];
    bn_muln_low(tmp, c2, k, 4);
    b2[0] = tmp[4];
    b2[1] = tmp[5];

	//Step 4: Compute parity corrections p1,p2 ahead of time
	//v2[0] is always even, so doesn't matter for parity.
	//b1*v1[0] is odd iff b1 is odd.
	p1 = (k[0] + b1[0] + 1) & 0x1;
	bn_add1_low(b1, b1, p1, 2);
	//v1[1] is always even, v2[1] odd.
	p2 = (b2[0] + 1) & 0x1;
    bn_add1_low(b2, b2, p2, 2);
	// (u1, u2) = (v1, v2) as alpha1 % 4 != 0

	//Step 5: Compute k1
    bn_muln_low(tmp, b1, v1[0], 2); //b1'*v1[0]
    bn_muln_low(tmp+4, b2, v2[0], 2); //b2'*v2[0]
    bn_addn_low(tmp+4, tmp, tmp+4, 4);
    *s1 = bn_subn_low(tmp, k, tmp+4, 4);
    //ADDACC_256(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7]); //b1'*v1[0] + b2'*v2[0]
	//SUBACC_256(tmp[4], tmp[5], tmp[6], tmp[7], k[0], k[1], k[2], k[3]); //k - (b1'*v1[0] + b2'*v2[0])
	//Take two's complement if needed.
	k1[0] = tmp[0] ^ (-*s1);
	k1[1] = tmp[1] ^ (-*s1);
	bn_add1_low(k1, k1, *s1, 2);

	//Step 6: Compute k2
	bn_muln_low(tmp, b1, v1[1], 2); //b1'*v1[1]
	bn_muln_low(tmp+4, b2, v2[1], 2); //-b2'*v2[1]
	//b1'*v1[1] > 0 and b2'*v2[1] < 0.
    *s2 = bn_subn_low(tmp, tmp+4, tmp, 4);
	//SUBACC_256(tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp2[0], tmp2[1], tmp2[2], tmp2[3]); //-b1'*v1[1] -b2'*v2[1]
	//Take two's complement if needed.
	k2[0] = tmp[0] ^ (-*s2);
	k2[1] = tmp[1] ^ (-*s2);
    bn_add1_low(k2, k2, *s2, 2);
}

/* 5-NAF pre-computation */
void smu_pre_4nf(__m128i *ppx0, __m128i *ppx1,
		__m128i *ppl0, __m128i *ppl1,
		__m128i *ppz0, __m128i *ppz1,
		__m128i px0, __m128i px1, __m128i pl0, __m128i pl1) {
	/* var */
	__m128i ONE, ml00;
	__m128i zin0[4], zin1[4];
	__m128i tmp0[1], tmp1[1];

	/* init */
	ONE = _mm_set_epi64x(0x0, 0x1);
	ml00 = _mm_xor_si128(pl0, ONE);

	/* pre-computation */
	/* P1 */
	ppx0[0] = px0;
	ppx1[0] = px1;
	ppl0[0] = pl0;
	ppl1[0] = pl1;

	/* P3 */
	eca_add_dbl(&ppx0[1], &ppx1[1], &ppl0[1], &ppl1[1], &ppz0[1], &ppz1[1],
			px0, px1, pl0, pl1, px0, px1, pl0, pl1);

	/* P5 and P7 */
	eca_add_sub_dbl(&ppx0[2], &ppx1[2], &ppl0[2], &ppl1[2], &ppz0[2], &ppz1[2],
			&ppx0[3], &ppx1[3], &ppl0[3], &ppl1[3], &ppz0[3], &ppz1[3],
			ppx0[1], ppx1[1], ppl0[1], ppl1[1], ppz0[1], ppz1[1],
			px0, px1, pl0, pl1);

	/* inversion: montgomery`s trick */
	/* part I */
	low_mul(&tmp0[0], &tmp1[0], ppz0[1], ppz1[1], ppz0[2], ppz1[2]);
	low_mul(&zin0[0], &zin1[0], tmp0[0], tmp1[0], ppz0[3], ppz1[3]);

	/* part II */
	low_inv(&zin0[0], &zin1[0], zin0[0], zin1[0]);

	/* part III */
	low_mul(&zin0[3], &zin1[3], zin0[0], zin1[0], tmp0[0], tmp1[0]);
	low_mul(&zin0[0], &zin1[0], zin0[0], zin1[0], ppz0[3], ppz1[3]);
	low_mul(&zin0[1], &zin1[1], zin0[0], zin1[0], ppz0[2], ppz1[2]);
	low_mul(&zin0[2], &zin1[2], zin0[0], zin1[0], ppz0[1], ppz1[1]);

	/* to affine */
	low_mul(&ppx0[1], &ppx1[1], ppx0[1], ppx1[1], zin0[1], zin1[1]);
	low_mul(&ppl0[1], &ppl1[1], ppl0[1], ppl1[1], zin0[1], zin1[1]);
	low_mul(&ppx0[2], &ppx1[2], ppx0[2], ppx1[2], zin0[2], zin1[2]);
	low_mul(&ppl0[2], &ppl1[2], ppl0[2], ppl1[2], zin0[2], zin1[2]);
	low_mul(&ppx0[3], &ppx1[3], ppx0[3], ppx1[3], zin0[3], zin1[3]);
	low_mul(&ppl0[3], &ppl1[3], ppl0[3], ppl1[3], zin0[3], zin1[3]);

	/* end */
	return;
}

/* 5-NAF pre-computation */
void smu_pre_5nf(__m128i *ppx0, __m128i *ppx1,
		__m128i *ppl0, __m128i *ppl1,
		__m128i *ppz0, __m128i *ppz1,
		__m128i px0, __m128i px1, __m128i pl0, __m128i pl1) {
	/* var */
	__m128i ONE, ml00;
	__m128i zin0[8], zin1[8];
	__m128i tmp0[5], tmp1[5];

	/* init */
	ONE = _mm_set_epi64x(0x0, 0x1);
	ml00 = _mm_xor_si128(pl0, ONE);

	/* pre-computation */
	/* P1 */
	ppx0[0] = px0;
	ppx1[0] = px1;
	ppl0[0] = pl0;
	ppl1[0] = pl1;

	/* P3 */
	eca_add_dbl(&ppx0[1], &ppx1[1], &ppl0[1], &ppl1[1], &ppz0[1], &ppz1[1],
			px0, px1, pl0, pl1, px0, px1, pl0, pl1);

	/* P5 and P7 */
	eca_add_sub_dbl(&ppx0[2], &ppx1[2], &ppl0[2], &ppl1[2], &ppz0[2], &ppz1[2],
			&ppx0[3], &ppx1[3], &ppl0[3], &ppl1[3], &ppz0[3], &ppz1[3],
			ppx0[1], ppx1[1], ppl0[1], ppl1[1], ppz0[1], ppz1[1],
			px0, px1, pl0, pl1);

	/* P9 and P11 */
	eca_add_sub_dbl(&ppx0[4], &ppx1[4], &ppl0[4], &ppl1[4], &ppz0[4], &ppz1[4],
			&ppx0[5], &ppx1[5], &ppl0[5], &ppl1[5], &ppz0[5], &ppz1[5],
			ppx0[2], ppx1[2], ppl0[2], ppl1[2], ppz0[2], ppz1[2],
			px0, px1, pl0, pl1);

	/* P9 and P11 */
	eca_add_sub_dbl(&ppx0[6], &ppx1[6], &ppl0[6], &ppl1[6], &ppz0[6], &ppz1[6],
			&ppx0[7], &ppx1[7], &ppl0[7], &ppl1[7], &ppz0[7], &ppz1[7],
			ppx0[3], ppx1[3], ppl0[3], ppl1[3], ppz0[3], ppz1[3],
			px0, px1, pl0, pl1);

	/* inversion: montgomery`s trick */
	/* part I */
	low_mul(&tmp0[0], &tmp1[0], ppz0[1], ppz1[1], ppz0[2], ppz1[2]);
	low_mul(&tmp0[1], &tmp1[1], tmp0[0], tmp1[0], ppz0[3], ppz1[3]);
	low_mul(&tmp0[2], &tmp1[2], tmp0[1], tmp1[1], ppz0[4], ppz1[4]);
	low_mul(&tmp0[3], &tmp1[3], tmp0[2], tmp1[2], ppz0[5], ppz1[5]);
	low_mul(&tmp0[4], &tmp1[4], tmp0[3], tmp1[3], ppz0[6], ppz1[6]);
	low_mul(&zin0[0], &zin1[0], tmp0[4], tmp1[4], ppz0[7], ppz1[7]);

	/* part II */
	low_inv_var(&zin0[0], &zin1[0], zin0[0], zin1[0]);

	/* part III */
	low_mul(&zin0[7], &zin1[7], zin0[0], zin1[0], tmp0[4], tmp1[4]);
	low_mul(&zin0[0], &zin1[0], zin0[0], zin1[0], ppz0[7], ppz1[7]);
	low_mul(&zin0[6], &zin1[6], zin0[0], zin1[0], tmp0[3], tmp1[3]);
	low_mul(&zin0[0], &zin1[0], zin0[0], zin1[0], ppz0[6], ppz1[6]);
	low_mul(&zin0[5], &zin1[5], zin0[0], zin1[0], tmp0[2], tmp1[2]);
	low_mul(&zin0[0], &zin1[0], zin0[0], zin1[0], ppz0[5], ppz1[5]);
	low_mul(&zin0[4], &zin1[4], zin0[0], zin1[0], tmp0[1], tmp1[1]);
	low_mul(&zin0[0], &zin1[0], zin0[0], zin1[0], ppz0[4], ppz1[4]);
	low_mul(&zin0[3], &zin1[3], zin0[0], zin1[0], tmp0[0], tmp1[0]);
	low_mul(&zin0[0], &zin1[0], zin0[0], zin1[0], ppz0[3], ppz1[3]);
	low_mul(&zin0[1], &zin1[1], zin0[0], zin1[0], ppz0[2], ppz1[2]);
	low_mul(&zin0[2], &zin1[2], zin0[0], zin1[0], ppz0[1], ppz1[1]);

	/* to affine */
	low_mul(&ppx0[1], &ppx1[1], ppx0[1], ppx1[1], zin0[1], zin1[1]);
	low_mul(&ppl0[1], &ppl1[1], ppl0[1], ppl1[1], zin0[1], zin1[1]);
	low_mul(&ppx0[2], &ppx1[2], ppx0[2], ppx1[2], zin0[2], zin1[2]);
	low_mul(&ppl0[2], &ppl1[2], ppl0[2], ppl1[2], zin0[2], zin1[2]);
	low_mul(&ppx0[3], &ppx1[3], ppx0[3], ppx1[3], zin0[3], zin1[3]);
	low_mul(&ppl0[3], &ppl1[3], ppl0[3], ppl1[3], zin0[3], zin1[3]);
	low_mul(&ppx0[4], &ppx1[4], ppx0[4], ppx1[4], zin0[4], zin1[4]);
	low_mul(&ppl0[4], &ppl1[4], ppl0[4], ppl1[4], zin0[4], zin1[4]);
	low_mul(&ppx0[5], &ppx1[5], ppx0[5], ppx1[5], zin0[5], zin1[5]);
	low_mul(&ppl0[5], &ppl1[5], ppl0[5], ppl1[5], zin0[5], zin1[5]);
	low_mul(&ppx0[6], &ppx1[6], ppx0[6], ppx1[6], zin0[6], zin1[6]);
	low_mul(&ppl0[6], &ppl1[6], ppl0[6], ppl1[6], zin0[6], zin1[6]);
	low_mul(&ppx0[7], &ppx1[7], ppx0[7], ppx1[7], zin0[7], zin1[7]);
	low_mul(&ppl0[7], &ppl1[7], ppl0[7], ppl1[7], zin0[7], zin1[7]);

	/* end */
	return;
}

/* retrieve digit sign (sig) and absolute value (abs) */
#define smu_get_flg(vec,idx,msk,abs,sig)\
    msk = vec[idx] >> 7;\
    abs = vec[idx];\
    sig = abs >> 63;\
    abs = ((abs ^ msk) + sig) >> 1;

/* GLS endomorphism */
#define smu_psi_end(ox00,ox01,ol00,ol01,ix00,ix01,il00,il01,ONE)\
    ox00 = _mm_srli_si128(ix00, 8);\
    ox00 = _mm_xor_si128(ox00, ix00);\
    ox01 = _mm_srli_si128(ix01, 8);\
    ox01 = _mm_xor_si128(ox01, ix01);\
    ol00 = _mm_srli_si128(il00, 8);\
    ol00 = _mm_xor_si128(ol00, il00);\
    ol01 = _mm_srli_si128(il01, 8);\
    ol01 = _mm_xor_si128(ol01, il01);\
    ol00 = _mm_xor_si128(ol00, ONE);

/* linar pass algorithm */
#define smu_lps(dst0, dst1, msk0, msk1, src, size)\
    dst0 = _mm_setzero_si128();\
    dst1 = _mm_setzero_si128();\
    for (int i = 0; i < size; i++) {\
        dst0 = _mm_xor_si128(dst0, _mm_and_si128(src[i], msk0[i]));\
        dst1 = _mm_xor_si128(dst1, _mm_and_si128(src[i], msk1[i]));\
    }

/* protected 5-NAF double-and-add left-to-right scalar multiplication */
void smu_5nf_dna_ltr(__m128i *qx0, __m128i *qx1, __m128i *ql0, __m128i *ql1,
		__m128i px0, __m128i px1, __m128i pl0, __m128i pl1, uint64_t *k) {
	/* var */
	int8_t dg0[32], dg1[32];
	int i, j;
	uint64_t k0[2], k1[2], sig0, sig1, abs0, abs1, msk, k0neg, k1neg;
	__m128i ppx0[8], ppl0[8], ppz0[8], ppx1[8], ppl1[8], ppz1[8];
	__m128i a0x0, a0x1, a0l0, a0l1, a1x0, a1x1, a1l0, a1l1;
	__m128i msk0[8], msk1[8], cmp[8], dig0, dig1, ssgn;
	__m128i qz0, qz1, e1x0, e1x1, e1l0, e1l1;

	/* init */
	__m128i ONE = _mm_set_epi64x(0x1, 0x1);
	cmp[0] = _mm_setzero_si128();
	for (j = 1; j < 8; j++) {
		cmp[j] = _mm_add_epi64(cmp[j - 1], ONE);
	}

	/* regular recoding */
	gls_recoding(k, k0, k1, &k0neg, &k1neg);
	smu_reg_rec(dg0, k0, 5);
	smu_reg_rec(dg1, k1, 5);

	/* pre computation */
	smu_pre_5nf(ppx0, ppx1, ppl0, ppl1, ppz0, ppz1, px0, px1, pl0, pl1);

	/* first iteration */
	/* digit */
	smu_get_flg(dg0, 31, msk, abs0, sig0);
	smu_get_flg(dg1, 31, msk, abs1, sig1);

    /* linear pass */
    dig0 = _mm_set_epi64x(abs0, abs0);
    dig1 = _mm_set_epi64x(abs1, abs1);
    for (j = 0; j < 8; j++) {
        msk0[j] = _mm_cmpeq_epi64(cmp[j], dig0);
        msk1[j] = _mm_cmpeq_epi64(cmp[j], dig1);
    }
    smu_lps(a0x0, a1x0, msk0, msk1, ppx0, 8);
    smu_lps(a0x1, a1x1, msk0, msk1, ppx1, 8);
    smu_lps(a0l0, a1l0, msk0, msk1, ppl0, 8);
    smu_lps(a0l1, a1l1, msk0, msk1, ppl1, 8);

    ssgn = _mm_set_epi64x(0x0, sig0 ^ k0neg);
    *qx0 = a0x0;
    *qx1 = a0x1;
	*ql0 = _mm_xor_si128(a0l0, ssgn);
	*ql1 = a0l1;

	/* add k1 digit */
	smu_psi_end(e1x0, e1x1, e1l0, e1l1, a1x0, a1x1, a1l0, a1l1, ONE);
	ssgn = _mm_set_epi64x(0x0, sig1 ^ k1neg);
    e1l0 = _mm_xor_si128(e1l0, ssgn);
	eca_add_mma(qx0, qx1, ql0, ql1, &qz0, &qz1, *qx0, *qx1, *ql0, *ql1, e1x0, e1x1, e1l0, e1l1);

	/* main loop */
	for (i = 30; i >= 0; i--) {
		/* point doubling */
		eca_dbl_ful(qx0, qx1, ql0, ql1, &qz0, &qz1,
				*qx0, *qx1, *ql0, *ql1, qz0, qz1);

		eca_dbl_ful(qx0, qx1, ql0, ql1, &qz0, &qz1,
				*qx0, *qx1, *ql0, *ql1, qz0, qz1);

		eca_dbl_ful(qx0, qx1, ql0, ql1, &qz0, &qz1,
				*qx0, *qx1, *ql0, *ql1, qz0, qz1);

		/* digit */
		smu_get_flg(dg0, i, msk, abs0, sig0);
		smu_get_flg(dg1, i, msk, abs1, sig1);

		/* linear pass */
		dig0 = _mm_set_epi64x(abs0, abs0);
		dig1 = _mm_set_epi64x(abs1, abs1);
		for (j = 0; j < 8; j++) {
			msk0[j] = _mm_cmpeq_epi64(cmp[j], dig0);
			msk1[j] = _mm_cmpeq_epi64(cmp[j], dig1);
		}
		smu_lps(a0x0, a1x0, msk0, msk1, ppx0, 8);
		smu_lps(a0x1, a1x1, msk0, msk1, ppx1, 8);
		smu_lps(a0l0, a1l0, msk0, msk1, ppl0, 8);
		smu_lps(a0l1, a1l1, msk0, msk1, ppl1, 8);

		/* add k0, k1 digits */
		ssgn = _mm_set_epi64x(0x0, sig0 ^ k0neg);
		a0l0 = _mm_xor_si128(a0l0, ssgn);

		smu_psi_end(e1x0, e1x1, e1l0, e1l1, a1x0, a1x1, a1l0, a1l1, ONE);
		ssgn = _mm_set_epi64x(0x0, sig1 ^ k1neg);
		e1l0 = _mm_xor_si128(e1l0, ssgn);

		eca_add_add_dbl(qx0, qx1, ql0, ql1, &qz0, &qz1, *qx0, *qx1, *ql0, *ql1,
                qz0, qz1, a0x0, a0x1, a0l0, a0l1, e1x0, e1x1, e1l0, e1l1);
	}

	/* to afffine */
	low_inv(&qz0, &qz1, qz0, qz1);
	low_mul(qx0, qx1, *qx0, *qx1, qz0, qz1);
	low_mul(ql0, ql1, *ql0, *ql1, qz0, qz1);

	/* final reduction */
	low_red_127_063_000(*qx0, *qx1, ONE);
	low_red_127_063_000(*ql0, *ql1, ONE);

	return;
}

/* protected 5-NAF double-and-add left-to-right scalar multiplication */
void smu_4nf_dna_ltr(__m128i *qx0, __m128i *qx1, __m128i *ql0, __m128i *ql1,
		__m128i px0, __m128i px1, __m128i pl0, __m128i pl1, uint64_t * k) {
	/* var */
	int8_t dg0[43], dg1[43];
	int i, j;
	uint64_t k0[2], k1[2], sig0, sig1, abs0, abs1, msk, k0neg, k1neg;
	__m128i ppx0[4], ppl0[4], ppz0[4], ppx1[4], ppl1[4], ppz1[4];
	__m128i a0x0, a0x1, a0l0, a0l1, a1x0, a1x1, a1l0, a1l1;
	__m128i msk0[4], msk1[4], cmp[4], dig0, dig1, ssgn;
	__m128i e1x0, e1x1, e1l0, e1l1, qz0, qz1;

	/* init */
	__m128i ONE = _mm_set_epi64x(0x1, 0x1);
	cmp[0] = _mm_setzero_si128();
	for (j = 1; j < 4; j++) {
		cmp[j] = _mm_add_epi64(cmp[j - 1], ONE);
	}

	/* regular recoding */
	gls_recoding(k, k0, k1, &k0neg, &k1neg);
	smu_reg_rec(dg0, k0, 4);
	smu_reg_rec(dg1, k1, 4);
    /* Negate endomorphism */
    k1neg ^= 1;

	/* pre computation */
	smu_pre_4nf(ppx0, ppx1, ppl0, ppl1, ppz0, ppz1, px0, px1, pl0, pl1);

    /* first iteration */
	/* digit */
	smu_get_flg(dg0, 42, msk, abs0, sig0);
	smu_get_flg(dg1, 42, msk, abs1, sig1);

    /* linear pass */
    dig0 = _mm_set_epi64x(abs0, abs0);
    dig1 = _mm_set_epi64x(abs1, abs1);
    for (j = 0; j < 4; j++) {
        msk0[j] = _mm_cmpeq_epi64(cmp[j], dig0);
        msk1[j] = _mm_cmpeq_epi64(cmp[j], dig1);
    }
    smu_lps(a0x0, a1x0, msk0, msk1, ppx0, 4);
    smu_lps(a0x1, a1x1, msk0, msk1, ppx1, 4);
    smu_lps(a0l0, a1l0, msk0, msk1, ppl0, 4);
    smu_lps(a0l1, a1l1, msk0, msk1, ppl1, 4);

    ssgn = _mm_set_epi64x(0x0, sig0 ^ k0neg);
    *qx0 = a0x0;
    *qx1 = a0x1;
	*ql0 = _mm_xor_si128(a0l0, ssgn);
	*ql1 = a0l1;

	/* add k1 digit */
	smu_psi_end(e1x0, e1x1, e1l0, e1l1, a1x0, a1x1, a1l0, a1l1, ONE);
	ssgn = _mm_set_epi64x(0x0, sig1 ^ k1neg);
    e1l0 = _mm_xor_si128(e1l0, ssgn);
	eca_add_mma(qx0, qx1, ql0, ql1, &qz0, &qz1, *qx0, *qx1, *ql0, *ql1, e1x0, e1x1, e1l0, e1l1);

	/* main loop */
	for (i = 41; i >= 0; i--) {
		/* point doubling */
		eca_dbl_ful(qx0, qx1, ql0, ql1, &qz0, &qz1,
				*qx0, *qx1, *ql0, *ql1, qz0, qz1);

		eca_dbl_ful(qx0, qx1, ql0, ql1, &qz0, &qz1,
				*qx0, *qx1, *ql0, *ql1, qz0, qz1);

		/* digit */
		smu_get_flg(dg0, i, msk, abs0, sig0);
		smu_get_flg(dg1, i, msk, abs1, sig1);

		/* linear pass */
		dig0 = _mm_set_epi64x(abs0, abs0);
		dig1 = _mm_set_epi64x(abs1, abs1);
		for (j = 0; j < 4; j++) {
			msk0[j] = _mm_cmpeq_epi64(cmp[j], dig0);
			msk1[j] = _mm_cmpeq_epi64(cmp[j], dig1);
		}
		smu_lps(a0x0, a1x0, msk0, msk1, ppx0, 4);
		smu_lps(a0x1, a1x1, msk0, msk1, ppx1, 4);
		smu_lps(a0l0, a1l0, msk0, msk1, ppl0, 4);
		smu_lps(a0l1, a1l1, msk0, msk1, ppl1, 4);

		/* add k0, k1 digits */
		ssgn = _mm_set_epi64x(0x0, sig0 ^ k0neg);
		a0l0 = _mm_xor_si128(a0l0, ssgn);
		smu_psi_end(e1x0, e1x1, e1l0, e1l1, a1x0, a1x1, a1l0, a1l1, ONE);
		ssgn = _mm_set_epi64x(0x0, sig1 ^ k1neg);
		e1l0 = _mm_xor_si128(e1l0, ssgn);
		eca_add_add_dbl(qx0, qx1, ql0, ql1, &qz0, &qz1,
				*qx0, *qx1, *ql0, *ql1, qz0, qz1,
				a0x0, a0x1, a0l0, a0l1, e1x0, e1x1, e1l0, e1l1);
	}

	/* to afffine */
	low_inv(&qz0, &qz1, qz0, qz1);
	low_mul(qx0, qx1, *qx0, *qx1, qz0, qz1);
	low_mul(ql0, ql1, *ql0, *ql1, qz0, qz1);

	/* final reduction */
	low_red_127_063_000(*qx0, *qx1, ONE);
	low_red_127_063_000(*ql0, *ql1, ONE);

	return;
}
