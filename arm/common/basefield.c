#include <stdio.h>
#include <stdlib.h>

#include "basefield.h"
#include "utils.h"

/**
 * basefield.c
 *
 * Implementation of arithmetic for the basefield F_{2^127}.
 * The field elements have a little endian representation,
 * meaning the least significant word has the lowest index.
 */

void bf_print_expr(poly64x2_t p) {
	poly64_t c;
	int wasFirst = 1;
	for(int i=1; i>=0; i--) {
		// 2^63 = the value of the leftmost bit in a word
		c = pow2to63;
		for(int j = 63; j>=0; j--) {
			poly64_t polybitcopy = c & p[i];
			if(polybitcopy == c) {
				if(!wasFirst) {
					printf("+");
				}
				wasFirst = 0;

				if(i == 0 && j == 0) {
					printf("1");
				} else {
					printf("z^%d", (i)*64 +j);
				}
			}
			c /= 2;
		}
	}

	if(wasFirst) {
		printf("0");
	}
}

void bf_print_expr_nl(poly64x2_t p) {
	bf_print_expr(p);
	printf("\n");
}

void bf_print_unred_expr(poly64x2x2_t p) {
	poly64_t c;
	int wasFirst = 1;
	for(int i=1; i>=0; i--) {
		for(int j=1; j>=0; j--) {
			// 2^63 = the value of the leftmost bit in a word
			c = pow2to63;
			for(int k = 63; k>=0; k--) {
				poly64_t polybitcopy = c & p.val[i][j];
				if(polybitcopy == c) {
					if(!wasFirst) {
						printf("+");
					}
					wasFirst = 0;

					if(i == 0 && j == 0 && k== 0) {
						printf("1");
					} else {
						printf("z^%d", i*128 +j*64 + k);
					}
				}
				c /= 2;
			}
		}
	}

	if(wasFirst) {
		printf("0");
	}
}

void bf_print_unred_expr_nl(poly64x2x2_t p) {
	bf_print_unred_expr(p);
	printf("\n");
}

void bf_print_hex(poly64x2_t p) {
	poly64x2_t red = bf_red_from_lazy(p);
	printf("%016lx||%016lx", red[1], red[0]);
}

void bf_print_hex_nl(poly64x2_t p) {
	bf_print_hex(p);
	printf("\n");
}

void bf_print_unred_hex(poly64x2x2_t p) {
	printf("%016lx||%016lx||%016lx||%016lx", p.val[1][1], p.val[1][0], p.val[0][1], p.val[0][0]);
}

void bf_print_unred_hex_nl(poly64x2x2_t p) {
	bf_print_unred_hex(p);
	printf("\n");
}

poly64x2_t bf_rand_elem() {
	// 2^63-1 = 01111111...
	long c = pow2to63-1;

	poly64_t p0 = rand_uint64();
	poly64_t p1 = rand_uint64() & c;

	poly64x2_t p = {p0, p1};

	return p;
}

poly64x2x2_t bf_pmull32(poly64x2_t a, poly64x2_t b) {
	poly64x2_t t;
	poly64x2x2_t r;
	r.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], b[0]));
	r.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[1], b[1]));
	t[1] = (poly64_t) veor_u64((uint64x1_t) b[0], (uint64x1_t) b[1]);
	t[0] = (poly64_t) veor_u64((uint64x1_t) a[0], (uint64x1_t) a[1]);
	t = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t[1], t[0]));
	t = (poly64x2_t) veorq_u64((uint64x2_t) t, (uint64x2_t) r.val[0]);
	t = (poly64x2_t) veorq_u64((uint64x2_t) t, (uint64x2_t) r.val[1]);
	r.val[0][1] = (poly64_t) veor_u64((uint64x1_t) r.val[0][1], (uint64x1_t) t[0]);
	r.val[1][0] = (poly64_t) veor_u64((uint64x1_t) r.val[1][0], (uint64x1_t) t[1]);

	return r;
}

//First tried to find an arm bit-interleaving instruction, but later realised this works:
poly64x2x2_t bf_psquare_neon(poly64x2_t a) {
	poly64x2x2_t r;
	r.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], a[0]));
	r.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[1], a[1]));
	return r;
}

/* Alg 2.40 - Modular reduction (one bit at a time) */

int has_reduction_precomputed = 0;

poly64x2_t reduction_polynomials[64];

void bf_red_generic_precomp() {
	if(has_reduction_precomputed) {
		return;
	}
	has_reduction_precomputed = 1;
	reduction_polynomials[0][0] = pow2to63+1;
	reduction_polynomials[0][1] = 0;
	for(int i = 1; i < 64; i++) {
		reduction_polynomials[i][0] = reduction_polynomials[i-1][0] << 1;
		reduction_polynomials[i][1] = reduction_polynomials[i-1][1] << 1;
		if(i==1) {
			reduction_polynomials[i][1]++;
		}
	}
}

poly64x2_t bf_red_generic(poly64x2x2_t c) {
	/* Step 1 */
	bf_red_generic_precomp();

	/* Step 2 */
	uint64_t digit_val = pow2to63 / 8;
	int index = 1;
	for(int digit = 252; digit >= 128; digit--) {
		/* Step 2.1 */
		if((c.val[1][index] & digit_val) == digit_val) {
			int j = (digit - 127)/64; //Can be 0 or 1
			int k = (digit - 127) - 64*j;
			if(j == 1) {
				c.val[0][1] = (poly64_t) veor_u64((uint64x1_t) reduction_polynomials[k][0], (uint64x1_t) c.val[0][1]);
				c.val[1][0] = (poly64_t) veor_u64((uint64x1_t) reduction_polynomials[k][1], (uint64x1_t) c.val[1][0]);
			} else {
				c.val[0] = bf_add(reduction_polynomials[k], c.val[0]);
			}
		}

		if(digit % 64 == 0) {
			digit_val = pow2to63;
			index--;
		} else {
			digit_val /= 2;
		}
	}
	if((c.val[0][1] & pow2to63) == pow2to63) {
		c.val[0][0] = reduction_polynomials[0][0] ^ c.val[0][0];
		c.val[0][1] &= pow2to63 - 1;
	}

	return c.val[0];
}

//For future ref: vshlq_n_u64 doesn't carry from uint64x2_t[0] to [1]
poly64x2_t bf_red_formula(poly64x2x2_t c) {
	poly64x2_t result = {0, 0};
	uint64_t bit127set = (c.val[0][1] & pow2to63) == pow2to63;
	uint64_t bit191set = (c.val[1][0] & pow2to63) == pow2to63;

	poly64x2_t term0to126 = {c.val[0][0], c.val[0][1] & (pow2to63 - 1)};
	result = bf_add(term0to126, result);

	poly64x2_t term252to127_rshift127 = {(c.val[1][0] << 1) + bit127set, (c.val[1][1] << 1) + bit191set};
	result = bf_add(term252to127_rshift127, result);

	poly64x2_t term190to127_rshift64 = {pow2to63*bit127set, c.val[1][0] & (pow2to63 - 1)};
	result = bf_add(term190to127_rshift64, result);

	poly64x2_t term252to191_rshift128 = {pow2to63*bit191set, c.val[1][1]};
	result = bf_add(term252to191_rshift128, result);

	poly64x2_t term252to191_rshift191 = {(c.val[1][1] << 1) + bit191set, 0};
	result = bf_add(term252to191_rshift191, result);

	return result;
}

poly64x2_t bf_red_neonv2(poly64x2x2_t c) {
	poly64x2_t a, t0, t1;
	
	t0 = vextq_p64(c.val[0], c.val[1], 1);
	t0[0] ^= c.val[1][0];
	t1 = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[1], 1);
	a[0] = (t0[0] << 1);
	a[0] >>= 1;
	a = vextq_p64(t1, a, 1);
	a = (poly64x2_t) veorq_u64((uint64x2_t) a, (uint64x2_t) t1);
	a[0] ^= c.val[0][0];
	t0 = (poly64x2_t) vshrq_n_u64((uint64x2_t) t0, 63);
	a = (poly64x2_t) veorq_u64((uint64x2_t) a, (uint64x2_t) t0);
	t0[0] <<= 63;
	t0 = vextq_p64(t0, t0, 1);
	t0 = (poly64x2_t) vzip2q_u64((uint64x2_t) t0, (uint64x2_t) c.val[1]);
	a = (poly64x2_t) veorq_u64((uint64x2_t) a, (uint64x2_t) t0);
	return a;
}

//Carry shift is not needed here!! 191 bit is always 0.
poly64x2_t bf_red_psquare_formula(poly64x2x2_t c) {
	poly64x2_t result = {0, 0};

	result = bf_add(c.val[0], result);

	poly64x2_t term255to128_rshift127 = {c.val[1][0] << 1, c.val[1][1] << 1};
	result = bf_add(term255to128_rshift127, result);

	poly64x2_t term191to128_rshift64 = {0, c.val[1][0]};
	result = bf_add(term191to128_rshift64, result);

	poly64x2_t term255to192_rshift128 = {0, c.val[1][1]};
	result = bf_add(term255to192_rshift128, result);

	poly64x2_t term255to192_rshift191 = {c.val[1][1] << 1, 0};
	result = bf_add(term255to192_rshift191, result);

	return result;
}

poly64x2_t bf_red_lazy_formula(poly64x2x2_t c) {
	poly64x2_t result = {0,0};
	
	result = bf_add(c.val[0], result);

	poly64x2_t term255to128_rshift127 = {c.val[1][0] << 1, c.val[1][1] << 1};
	result = bf_add(term255to128_rshift127, result);
	
	poly64x2_t bit191 = {0, c.val[1][0] >> 63};
	result = bf_add(bit191, result);

	poly64x2_t term191to128_rshift64 = {0, c.val[1][0]};
	result = bf_add(term191to128_rshift64, result);

	poly64x2_t term255to192_rshift128 = {0, c.val[1][1]};
	result = bf_add(term255to192_rshift128, result);

	poly64x2_t term255to192_rshift191 = {c.val[1][1] << 1, 0};
	result = bf_add(term255to192_rshift191, result);

	return result;
}

poly64x2_t bf_inv(poly64x2_t a) {
	poly64x2_t x_10 = bf_red_psquare(bf_psquare(a)); // 2
	poly64x2_t x_11 = bf_red(bf_pmull(a, x_10)); //1 + 2 = 3
	poly64x2_t x_110 = bf_red_psquare(bf_psquare(x_11)); //3*2 = 6
	poly64x2_t x_111 = bf_red(bf_pmull(a, x_110)); //1 + 6 = 7
	poly64x2_t x_111000 = bf_red_psquare(bf_psquare(bf_red_psquare(bf_psquare(bf_red_psquare(bf_psquare(x_111))))));
	//7*2^3 = 56
	poly64x2_t x_111111 = bf_red(bf_pmull(x_111, x_111000)); //56 + 7 = 2^6 - 1
	poly64x2_t x_x12 = bf_red(bf_pmull(bf_multisquare_loop(x_111111, 6), x_111111));
	//(2^6 -1)*2^6 + 2^6 - 1 = 2^12 - 1
	poly64x2_t x_x24 = bf_red(bf_pmull(bf_multisquare_loop(x_x12, 12), x_x12));
	//(2^12-1)*2^12 + 2^12 - 1 = 2^24 - 1
	poly64x2_t x_i34 = bf_multisquare_loop(x_x24, 6); //(2^24-1)*2^6 = 2^30 -2^6
	poly64x2_t x_x30 = bf_red(bf_pmull(x_111111, x_i34)); //(2^30-2^6) + 2^6 - 1 = 2^30 - 1
	poly64x2_t x_x48 = bf_red(bf_pmull(bf_multisquare_loop(x_i34, 18), x_x24));
	//(2^30 - 2^6)*2^18 + 2^24 - 1 = 2^48 - 1
	poly64x2_t x_x96 = bf_red(bf_pmull(bf_multisquare_loop(x_x48, 48), x_x48));
	//(2^48-1)*2^48 +2^48 - 1 = 2^96 -1
	poly64x2_t x_end = bf_red(bf_pmull(bf_multisquare_loop(x_x96, 30), x_x30));
	//(2^96-1)*2^30 + 2^30 - 1 = 2^126 - 1
	return bf_red_psquare(bf_psquare(x_end));
}

poly64x2x2_t bf_pmull64(poly64x2_t a, poly64x2_t b) {
	poly64x2_t t0, t1;
	poly64x2_t z = {0,0};
	poly64x2x2_t r;

	r.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], b[0]));
	r.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a, b));
	t0 = vextq_p64(b, b, 1);
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], t0[0]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a, t0));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	r.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) r.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	r.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) r.val[1], (uint64x2_t) t1);

	return r;

}

//veor3q not defined?
poly64x2_t bf_red_psquare_neon(poly64x2x2_t c) {
	c.val[1][0] = (poly64_t) veor_u64((uint64x1_t) c.val[1][0], (uint64x1_t) c.val[1][1]);
	c.val[0][1] = (poly64_t) veor_u64((uint64x1_t) c.val[0][1], (uint64x1_t) c.val[1][0]);
	c.val[1] = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[1], 1);
	return (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) c.val[1]);
}

poly64x2_t bf_red_psquare_neonv2(poly64x2x2_t c) {
	poly64x2_t t;
	t[0] = 0;
	t = vextq_p64(c.val[1], t, 1);
	c.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[1], (uint64x2_t) t);
	t = vextq_p64(t, c.val[1], 1);
	c.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) t);
	c.val[1] = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[1], 1);
	return (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) c.val[1]);
}

poly64x2_t bf_red_from_lazy(poly64x2_t a) {
	uint64_t t = a[1] & 0x8000000000000000;
	a[1] ^= t;
	t ^= (t >> 63);
	a[0] ^= t;
	return a;
} 

//Simple and slow, compute a^(-1) as a^((2^127)-2).
//Exploits the fact that 2^127 - 1 = 2^126 + 2^125 + ... + 2^1 + 1,
//hence a^-1 = a^(2^126)*...*a^(2^3)*a^(2^2)*a^2:
poly64x2_t bf_fermat_inv(poly64x2_t a) {
	poly64x2_t inv = {1, 0};
	poly64x2_t power = {a[0], a[1]};
	for(int i = 0; i < 126; i++) {
		power = bf_red_psquare(bf_psquare(power));
		inv = bf_red(bf_pmull(inv, power));
	}
	return inv;
}

//Addition chain for 126
//Addition chain means next term is sum of two previous terms.
// 1 -> 2 -> 3 -> 6 -> 12 -> 24 -> 30 -> 48 -> 96 -> 126
//Means need 9 multiplications & 126 squarings using Itoh & Tsujii alg
//In the end return (a^(2^126 -1))^2 = a^(2^127 -2) = a^-1 per Fermat
poly64x2_t bf_addchain_inv(poly64x2_t a) {
	poly64x2_t x_10 = bf_red_psquare(bf_psquare(a)); // 2
	poly64x2_t x_11 = bf_red_lazy(bf_pmull(a, x_10)); //1 + 2 = 3
	poly64x2_t x_110 = bf_red_psquare(bf_psquare(x_11)); //3*2 = 6
	poly64x2_t x_111 = bf_red_lazy(bf_pmull(a, x_110)); //1 + 6 = 7
	poly64x2_t x_111000 = bf_red_psquare(bf_psquare(bf_red_psquare(bf_psquare(bf_red_psquare(bf_psquare(x_111))))));
	//7*2^3 = 56
	poly64x2_t x_111111 = bf_red_lazy(bf_pmull(x_111, x_111000)); //56 + 7 = 2^6 - 1
	poly64x2_t x_x12 = bf_red_lazy(bf_pmull(bf_multisquare_loop(x_111111, 6), x_111111));
	//(2^6 -1)*2^6 + 2^6 - 1 = 2^12 - 1
	poly64x2_t x_x24 = bf_red_lazy(bf_pmull(bf_multisquare_loop(x_x12, 12), x_x12));
	//(2^12-1)*2^12 + 2^12 - 1 = 2^24 - 1
	poly64x2_t x_i34 = bf_multisquare_loop(x_x24, 6); //(2^24-1)*2^6 = 2^30 -2^6
	poly64x2_t x_x30 = bf_red_lazy(bf_pmull(x_111111, x_i34)); //(2^30-2^6) + 2^6 - 1 = 2^30 - 1
	poly64x2_t x_x48 = bf_red_lazy(bf_pmull(bf_multisquare_loop(x_i34, 18), x_x24));
	//(2^30 - 2^6)*2^18 + 2^24 - 1 = 2^48 - 1
	poly64x2_t x_x96 = bf_red_lazy(bf_pmull(bf_multisquare_loop(x_x48, 48), x_x48));
	//(2^48-1)*2^48 +2^48 - 1 = 2^96 -1
	poly64x2_t x_end = bf_red_lazy(bf_pmull(bf_multisquare_loop(x_x96, 30), x_x30));
	//(2^96-1)*2^30 + 2^30 - 1 = 2^126 - 1
	return bf_red_psquare(bf_psquare(x_end));
}

//Infinite loops with gcc???
//Algorithm 2.49
poly64x2_t bf_nonconst_inv(poly64x2_t a) {
	poly64x2_t f = {pow2to63 + 1, pow2to63};
	poly64x2_t u = a;
	poly64x2_t v = f;
	poly64x2_t g1 = {1, 0};
	poly64x2_t g2 = {0, 0};
	int uIsOne = (u[0] == 1 && u[1] == 0);
	int vIsOne = 0;
	//int i = 0;
	while(!uIsOne && !vIsOne) {
		//printf("i = %d\n", i);
		//i++;
		//printf("u ");
		//bf_print_hex_nl(u);
		//printf("v ");
		//bf_print_hex_nl(v);
		while(!(u[0] == 0 && u[1] == 0) && u[0] % 2 == 0) {
			//printf("u[0] = %lu divisible by z\n", u[0]);
			u = shift_right_carry(u);
			if(g1[0] % 2 == 0) {
				g1 = shift_right_carry(g1);
			} else {
				g1 = bf_add(g1, f);
				g1 = shift_right_carry(g1);
			}
		}
		while(!(v[0] == 0 && v[1] == 0) && v[0] % 2 == 0) {
			//printf("v[0] = %lu divisible by z\n", v[0]);
			v = shift_right_carry(v);
			if(g2[0] % 2 == 0) {
				g2 = shift_right_carry(g2);
			} else {
				g2 = bf_add(g2, f);
				g2 = shift_right_carry(g2);
			}
		}

		uint64_t greaterDeg = hasGreaterDeg(u, v);
		//printf("greater deg: %lu\n", greaterDeg);
		if(greaterDeg) {
			u = bf_add(u, v);
			g1 = bf_add(g1, g2);
		} else {
			v = bf_add(u,v);
			g2 = bf_add(g1, g2);
		}

		uIsOne = (u[0] == 1 && u[1] == 0);
		vIsOne = (v[0] == 1 && v[1] == 0);
	}

	if(uIsOne) return g1;
	return g2;
}

//Need multisquaring by 6, 12, 18, 30, 48, will use precomp approach
//from "2 is the fastest prime".
//For fixed k, need table of dim 32x16
poly64x2_t T_6[32][16];
poly64x2_t T_12[32][16];
poly64x2_t T_18[32][16];
poly64x2_t T_30[32][16];
poly64x2_t T_48[32][16];
uint64_t has_precomputed_inv_tables = 0;

void precomp_inv_table(uint64_t k) {
	poly64x2_t zpow_4 = {16, 0};
	poly64x2_t zpow_4j = {1, 0};
	poly64x2_t zpow_4jplus1 = {2, 0};
	poly64x2_t zpow_4jplus2 = {4, 0};
	poly64x2_t zpow_4jplus3 = {8, 0};
	for(int j = 0; j < 32; j++) {
		for(int i = 0; i < 16; i++) {
			poly64x2_t sum = {0, 0};
			if(i % 2 == 1) { //bit i0 is set
				sum = bf_add(sum, zpow_4j);
			}
			uint64_t l = i / 2;
			if(l % 2 == 1) { //bit i1 is set
				sum = bf_add(sum, zpow_4jplus1);
			}
			l /= 2;
			if(l % 2 == 1) { //bit i2 is set
				sum = bf_add(sum, zpow_4jplus2);
			}
			l /= 2;
			if(l % 2 == 1) { //bit i3 is set
				sum = bf_add(sum, zpow_4jplus3);
			}
			if(k == 6) {
				T_6[j][i] = bf_multisquare_loop(sum, 6);
			} else if(k == 12) {
				T_12[j][i] = bf_multisquare_loop(sum, 12);
			} else if(k == 18) {
				T_18[j][i] = bf_multisquare_loop(sum, 18);
			} else if(k == 30) {
				T_30[j][i] = bf_multisquare_loop(sum, 30);
			} else if(k == 48) {
				T_48[j][i] = bf_multisquare_loop(sum, 48);
			}
		}
		zpow_4j = bf_red(bf_pmull(zpow_4j, zpow_4));
		zpow_4jplus1 = bf_red(bf_pmull(zpow_4jplus1, zpow_4));
		zpow_4jplus2 = bf_red(bf_pmull(zpow_4jplus2, zpow_4));
		zpow_4jplus3 = bf_red(bf_pmull(zpow_4jplus3, zpow_4));
	}
}

void precomp_inv_tables() {
	if(has_precomputed_inv_tables) {
		return;
	}

	// precomp_inv_table(6);
	// precomp_inv_table(12);
	//precomp_inv_table(18);
	precomp_inv_table(30);
	precomp_inv_table(48);

	has_precomputed_inv_tables = 1;
}

// void free_inv_tables() {
// 	free(T_6);
// 	free(T_12);
// 	free(T_18);
// 	free(T_30);
// 	free(T_48);
// }

poly64x2_t bf_multisquare_lookup_6(poly64x2_t a) {
	poly64x2_t res = (poly64x2_t) veorq_u64((uint64x2_t) T_6[0][a[0] & 15], (uint64x2_t) T_6[16][a[1] & 15]);
	for(int j = 1; j < 16; j++) {
		a = (poly64x2_t) vshrq_n_u64((uint64x2_t) a, 4);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_6[j][a[0] & 15]);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_6[j+16][a[1] & 15]);
	}
	return res;
}

poly64x2_t bf_multisquare_lookup_12(poly64x2_t a) {
	poly64x2_t res = (poly64x2_t) veorq_u64((uint64x2_t) T_12[0][a[0] & 15], (uint64x2_t) T_12[16][a[1] & 15]);
	for(int j = 1; j < 16; j++) {
		a = (poly64x2_t) vshrq_n_u64((uint64x2_t) a, 4);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_12[j][a[0] & 15]);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_12[j+16][a[1] & 15]);
	}
	return res;
}

poly64x2_t bf_multisquare_lookup_18(poly64x2_t a) {
	poly64x2_t res = (poly64x2_t) veorq_u64((uint64x2_t) T_18[0][a[0] & 15], (uint64x2_t) T_18[16][a[1] & 15]);
	for(int j = 1; j < 16; j++) {
		a = (poly64x2_t) vshrq_n_u64((uint64x2_t) a, 4);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_18[j][a[0] & 15]);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_18[j+16][a[1] & 15]);
	}
	return res;
}

poly64x2_t bf_multisquare_lookup_30(poly64x2_t a) {
	poly64x2_t res = (poly64x2_t) veorq_u64((uint64x2_t) T_30[0][a[0] & 15], (uint64x2_t) T_30[16][a[1] & 15]);
	for(int j = 1; j < 16; j++) {
		a = (poly64x2_t) vshrq_n_u64((uint64x2_t) a, 4);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_30[j][a[0] & 15]);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_30[j+16][a[1] & 15]);
	}
	return res;
}

poly64x2_t bf_multisquare_lookup_48(poly64x2_t a) {
	poly64x2_t res = (poly64x2_t) veorq_u64((uint64x2_t) T_48[0][a[0] & 15], (uint64x2_t) T_48[16][a[1] & 15]);
	for(int j = 1; j < 16; j++) {
		a = (poly64x2_t) vshrq_n_u64((uint64x2_t) a, 4);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_48[j][a[0] & 15]);
		res = (poly64x2_t) veorq_u64((uint64x2_t) res, (uint64x2_t) T_48[j+16][a[1] & 15]);
	}
	return res;
}

poly64x2_t bf_addchain_lookup_inv(poly64x2_t a) {
	poly64x2_t t1 = bf_red_psquare(bf_psquare(a)); // 2
	t1 = bf_red(bf_pmull(a, t1)); //1 + 2 = 3
	t1 = bf_red_psquare(bf_psquare(t1)); //3*2 = 6
	t1 = bf_red(bf_pmull(a, t1)); //1 + 6 = 7
	a = bf_red_psquare(bf_psquare(bf_red_psquare(bf_psquare(bf_red_psquare(bf_psquare(t1))))));
	//7*2^3 = 56
	a = bf_red(bf_pmull(t1, a)); //56 + 7 = 2^6 - 1
	t1 = bf_red(bf_pmull(bf_multisquare_loop(a, 6), a));
	//(2^6 -1)*2^6 + 2^6 - 1 = 2^12 - 1
	t1 = bf_red(bf_pmull(bf_multisquare_loop(t1, 12), t1));
	//(2^12-1)*2^12 + 2^12 - 1 = 2^24 - 1
	poly64x2_t t2 = bf_multisquare_loop(t1, 6); //(2^24-1)*2^6 = 2^30 -2^6
	a = bf_red(bf_pmull(a, t2)); //(2^30-2^6) + 2^6 - 1 = 2^30 - 1
	t1 = bf_red(bf_pmull(bf_multisquare_loop(t2, 18), t1));
	//(2^30 - 2^6)*2^18 + 2^24 - 1 = 2^48 - 1
	t1 = bf_red(bf_pmull(bf_multisquare_lookup_48(t1), t1));
	//(2^48-1)*2^48 +2^48 - 1 = 2^96 -1
	a = bf_red(bf_pmull(bf_multisquare_lookup_30(t1), a));
	//(2^96-1)*2^30 + 2^30 - 1 = 2^126 - 1
	return bf_red_psquare(bf_psquare(a));
}

poly64x2_t bf_red_lazy1(poly64x2x2_t c) {
	poly64x2_t t1, t2;
	t1[0] = 0;
	t2[0] = c.val[1][0] >> 63;
	t1 = vextq_p64(c.val[1], t1, 1);
	c.val[1][0] ^= t1[0];
	t2[0] ^= c.val[1][0];
	t1 = vextq_p64(t1, t2, 1);
	c.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) t1);
	c.val[1] = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[1], 1);
	return (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) c.val[1]);
}
