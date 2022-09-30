/* inc */
#include <assert.h>
#include <stdio.h>
#include "lib.h"
#include "api.h"

#include "bn.h"
#include "ffa.h"
#include "eca.h"
#include "enc.h"
#include "smu.h"

/* Return the generator. */
int crypto_dh_generator(unsigned char *pk) {
	__m128i px0, px1, pl0, pl1;
	px0 = _mm_set_epi64x(0x9D1932CB5FA5B9BF, 0x5BE5F4EB93D8712A);
	px1 = _mm_set_epi64x(0x25F2F29FCBDEC78E, 0x47E70D2DCA8C7210);
	pl0 = _mm_set_epi64x(0x25BE90C01E0E9B06, 0x97FBBBBFEB3A8AB4);
	pl1 = _mm_set_epi64x(0x0B3834B048C217C1, 0x1A1764D658204447);
	ec_enc(pk, px0, px1, pl0, pl1);
	return 0;
}

/* Key pair generation. */
int crypto_dh_gls254prot_opt_keypair(unsigned char *pk, unsigned char *sk) {
	__m128i px0, px1, pl0, pl1;

	/* Mask top bits to reduce scalar modulo order. */
	sk[31] &= 0x3F;

	/* generator */
	px0 = _mm_set_epi64x(0x9D1932CB5FA5B9BF, 0x5BE5F4EB93D8712A);
	px1 = _mm_set_epi64x(0x25F2F29FCBDEC78E, 0x47E70D2DCA8C7210);
	pl0 = _mm_set_epi64x(0x25BE90C01E0E9B06, 0x97FBBBBFEB3A8AB4);
	pl1 = _mm_set_epi64x(0x0B3834B048C217C1, 0x1A1764D658204447);

	/* Scalar multiplication. */
	smu_5nf_dna_ltr(&px0, &px1, &pl0, &pl1, px0, px1, pl0, pl1, (uint64_t *)sk);

	/* Write the result. */
	ec_enc(pk, px0, px1, pl0, pl1);

	return 0;
}

/* Shared secret computation. */
int crypto_dh_gls254prot_1dw4(unsigned char *out, unsigned char *pk,
                unsigned char *sk) {
        __m128i px0, px1, pl0, pl1;

        sk[31] &= 0x3F;

        ec_dec(&px0, &px1, &pl0, &pl1, pk);

        smu_4nf_dna_ltr(&px0, &px1, &pl0, &pl1, px0, px1, pl0, pl1, (uint64_t *)sk);

        ec_enc(out, px0, px1, pl0, pl1);

        return 0;
}

/* Shared secret computation. */
int crypto_dh_gls254prot_1dw5(unsigned char *out, unsigned char *pk,
		unsigned char *sk) {
	__m128i px0, px1, pl0, pl1;

	sk[31] &= 0x3F;

	ec_dec(&px0, &px1, &pl0, &pl1, pk);

	smu_5nf_dna_ltr(&px0, &px1, &pl0, &pl1, px0, px1, pl0, pl1, (uint64_t *)sk);

	ec_enc(out, px0, px1, pl0, pl1);

	return 0;
}

/* Shared secret computation. */
int crypto_dh_gls254prot_2dw3(unsigned char *out, unsigned char *pk,
		unsigned char *sk) {
	__m128i px0, px1, pl0, pl1;

	sk[31] &= 0x3F;

	ec_dec(&px0, &px1, &pl0, &pl1, pk);

	smu_3nf_2d_ltr(&px0, &px1, &pl0, &pl1, px0, px1, pl0, pl1, (uint64_t *)sk);

	ec_enc(out, px0, px1, pl0, pl1);

	return 0;
}

/* Shared secret computation. */
int crypto_dh_gls254prot_2dw4(unsigned char *out, unsigned char *pk,
		unsigned char *sk) {
	__m128i px0, px1, pl0, pl1;

	sk[31] &= 0x3F;

	ec_dec(&px0, &px1, &pl0, &pl1, pk);

	smu_4nf_2d_ltr(&px0, &px1, &pl0, &pl1, px0, px1, pl0, pl1, (uint64_t *)sk);

	ec_enc(out, px0, px1, pl0, pl1);

	return 0;
}

int crypto_dh_gls254prot_hash(unsigned char *out, unsigned char *sk) {
	__m128i x0, x1, l0, l1, z0, z1, u0, u1;
	__m128i _x0, _x1, _l0, _l1;
	sk[31] &= 0x7F;
	u0 = _mm_loadu_si128((__m128i *) & sk[0]);
	u1 = _mm_loadu_si128((__m128i *) & sk[16]);
	ec_sw(&x0, &x1, &l0, &l1, u0, u1);
	ec_sw(&_x0, &_x1, &_l0, &_l1, u1, u0);
	eca_add_mma(&x0, &x1, &l0, &l1, &z0, &z1, x0, x1, l0, l1, _x0, _x1, _l0, _l1);
	eca_dbl_ful(&x0, &x1, &l0, &l1, &z0, &z1, x0, x1, l0, l1, z0, z1);

    low_inv(&z0, &z1, z0, z1);
    low_mul(&x0, &x1, x0, x1, z0, z1);
    low_mul(&l0, &l1, l0, l1, z0, z1);
    low_red_127_063_000(x0, x1, z0);
    low_red_127_063_000(l0, l1, z0);
    ec_enc(out, x0, x1, l0, l1);
	return 0;
}

#ifdef MAIN

#include <string.h>
#include "bench.h"
#include "bench.c"

static void ec_test() {
	uint8_t p[64], q[64];
	unsigned long long int u[4];
	__m128i x0, x1, l0, l1;

	u[0] = 0x3CBDE37CF43A8CF4;
	u[1] = 0x3F1A47DEDC1A1DAD;
	u[2] = 0x0;
	u[3] = 0x2000000000000000;

	crypto_dh_generator(q);
	crypto_dh_gls254prot_2dw3(p, q, (unsigned char*)u);
	assert(memcmp(p, q, 32) == 0);


	for (int i = 0; i < 4; i++) {
		__builtin_ia32_rdrand64_step(&u[i]);
	}

	crypto_dh_generator(q);
	crypto_dh_gls254prot_opt_keypair(p, (unsigned char*)u);
	crypto_dh_gls254prot_1dw5(q, q, (unsigned char*)u);
	assert(memcmp(p, q, 32) == 0);
	crypto_dh_generator(q);
	crypto_dh_gls254prot_2dw3(q, q, (unsigned char*)u);
	assert(memcmp(p, q, 32) == 0);
	crypto_dh_generator(q);
	crypto_dh_gls254prot_2dw4(q, q, (unsigned char*)u);
	assert(memcmp(p, q, 32) == 0);

	l0 = _mm_loadu_si128((__m128i *) u);
	l1 = _mm_loadu_si128((__m128i *) (u + 2));
	ec_sw(&x0, &x1, &l0, &l1, l0, l1);
	assert(ec_ok(x0, x1, l0, l1) == 1);
	/* Multiply by cofactor. */
	eca_dbl_aff(&x0, &x1, &l0, &l1, x0, x1, l0, l1);
	ec_enc(p, x0, x1, l0, l1);
	assert(ec_dec(&x0, &x1, &l0, &l1, p));

	crypto_dh_gls254prot_opt_keypair(p, (unsigned char*)u);
	assert(ec_dec(&x0, &x1, &l0, &l1, p));
	crypto_dh_gls254prot_2dw3(p, p, (unsigned char *)u);
	assert(ec_dec(&x0, &x1, &l0, &l1, p));
	ec_enc(q, x0, x1, l0, l1);
	assert(memcmp(p, q, 32) == 0);
	memset(p, 0, sizeof(p));
	assert(ec_dec(&x0, &x1, &l0, &l1, p) == 0);
	assert(ec_ok(x0, x1, l0, l1) == 0);
}

static void dh_test() {
	uint8_t pa[64], pb[64], k1[64], k2[64];
	unsigned long long int sa[4], sb[4];

	for (int i = 0; i < 4; i++) {
		__builtin_ia32_rdrand64_step(&sa[i]);
		__builtin_ia32_rdrand64_step(&sb[i]);
	}

	crypto_dh_gls254prot_opt_keypair(pa, (unsigned char *)sa);
	crypto_dh_gls254prot_opt_keypair(pb, (unsigned char *)sb);
	crypto_dh_gls254prot_2dw3(k1, pa, (unsigned char *)sb);
	crypto_dh_gls254prot_2dw3(k2, pb, (unsigned char *)sa);
	assert(memcmp(k1, k2, 32) == 0);
}

void bench() {
	uint8_t p[64], q[64];
	unsigned long long int u[4], l[4];
	__m128i x0, x1, l0, l1;

	printf("-- Low-level benchmarks:\n\n");

	BENCH_BEGIN("low_sqr_bas") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		x0 = _mm_loadu_si128((__m128i *) u);
		BENCH_ADD(low_sqr_bas(&x0, x0));
	} BENCH_END;

	BENCH_BEGIN("low_mul_bas") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		x0 = _mm_loadu_si128((__m128i *) u);
		x1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_mul_bas(&x0, x0, x1));
	} BENCH_END;

	BENCH_BEGIN("low_red_bas") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		l0 = _mm_loadu_si128((__m128i *) u);
		l1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_red_128_064_001_bas(x0, x1, l0, l1, x0, x1));
	} BENCH_END;

	BENCH_BEGIN("low_inv_bas") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		x0 = _mm_loadu_si128((__m128i *) u);
		BENCH_ADD(low_inv_bas(&x0, x0));
	} BENCH_END;

	BENCH_BEGIN("low_inv_tbl") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		x0 = _mm_loadu_si128((__m128i *) u);
		BENCH_ADD(low_inv_tbl(&x0, x0));
	} BENCH_END;

	BENCH_BEGIN("low_sqr") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		l0 = _mm_loadu_si128((__m128i *) u);
		l1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_sqr(&x0, &x1, x0, x1));
	} BENCH_END;

	BENCH_BEGIN("low_mul") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		x0 = _mm_loadu_si128((__m128i *) u);
		x1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_mul(&x0, &x1, x0, x1, l0, l1));
	} BENCH_END;

	BENCH_BEGIN("low_red") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		l0 = _mm_loadu_si128((__m128i *) u);
		l1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_red_128_064_001(x0, x1, l0, l1, x0, x1, x0, x1));
	} BENCH_END;

	BENCH_BEGIN("low_inv") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		l0 = _mm_loadu_si128((__m128i *) u);
		l1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_inv(&x0, &x1, l0, l1));
	} BENCH_END;

	BENCH_BEGIN("low_inv_var") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		x0 = _mm_loadu_si128((__m128i *) u);
		x1 = _mm_loadu_si128((__m128i *) (u + 2));
		BENCH_ADD(low_inv_var(&x0, &x1, x0, x1));
	} BENCH_END;

	BENCH_BEGIN("ec_sw") {
		BENCH_ADD(ec_sw(&x0, &x1, &l0, &l1, l0, l1));
	} BENCH_END;

	printf("-- Scalar multiplication benchmarks:\n\n");

	BENCH_BEGIN("crypto_dh_gls254prot_opt_keypair") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		BENCH_ADD(crypto_dh_gls254prot_opt_keypair(p, (unsigned char *)u));
	} BENCH_END;

	BENCH_BEGIN("crypto_dh_gls254prot_1dw4") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		BENCH_ADD(crypto_dh_gls254prot_1dw4(p, p, (unsigned char *)u));
	} BENCH_END;

	BENCH_BEGIN("crypto_dh_gls254prot_1dw5") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		BENCH_ADD(crypto_dh_gls254prot_1dw5(p, p, (unsigned char *)u));
	} BENCH_END;

	BENCH_BEGIN("crypto_dh_gls254prot_2dw3") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		BENCH_ADD(crypto_dh_gls254prot_2dw3(p, p, (unsigned char *)u));
	} BENCH_END;

	BENCH_BEGIN("crypto_dh_gls254prot_2dw4") {
		for (int i = 0; i < 4; i++) {
			__builtin_ia32_rdrand64_step(&u[i]);
		}
		BENCH_ADD(crypto_dh_gls254prot_2dw4(p, p, (unsigned char *)u));
	} BENCH_END;
}

int main(int argc, char const *argv[]) {
	ec_ell_pre();
	for (int i = 0; i < 10000; i++) {
		ec_test();
		dh_test();
	}

	printf("\nTests PASSED!\n\n");

	bench();
	return 0;
}

#endif
