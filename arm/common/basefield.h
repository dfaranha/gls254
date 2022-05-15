#include <arm_neon.h>
//#include <stdio.h>

#ifndef BASEFIELD_H
#define BASEFIELD_H

#define pow2to63 9223372036854775808U

void bf_print_expr(poly64x2_t p);

void bf_print_expr_nl(poly64x2_t p);

void bf_print_unred_expr(poly64x2x2_t p);

void bf_print_unred_expr_nl(poly64x2x2_t p);

void bf_print_hex(poly64x2_t p);

void bf_print_hex_nl(poly64x2_t p);

void bf_print_unred_hex(poly64x2x2_t p);

void bf_print_unred_hex_nl(poly64x2x2_t p);

static inline poly64x2_t bf_create_elem(uint64_t l, uint64_t h) {
	poly64x2_t a = {l, h};
	return a;
}

poly64x2_t bf_rand_elem();

static inline poly64x2_t bf_add(poly64x2_t a, poly64x2_t b) {
	return (poly64x2_t) veorq_u64((uint64x2_t) a, (uint64x2_t) b);
}

static inline poly64x2x2_t bf_pmull(poly64x2_t a, poly64x2_t b) {
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

static inline poly64x2x2_t bf_psquare(poly64x2_t a) {
	poly64x2x2_t r;
	r.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], a[0]));
	r.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[1], a[1]));
	return r;
}

static inline poly64x2_t bf_red(poly64x2x2_t c) {
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
	t0 = (poly64x2_t) vextq_u64((uint64x2_t) t0, (uint64x2_t) t0, 1);
	t0 = (poly64x2_t) vzip2q_u64((uint64x2_t) t0, (uint64x2_t) c.val[1]);
	a = (poly64x2_t) veorq_u64((uint64x2_t) a, (uint64x2_t) t0);
	return a;
}

static inline poly64x2_t bf_red_psquare(poly64x2x2_t c) {
	poly64x2_t t;
	t[0] = 0;
	t = vextq_p64(c.val[1], t, 1);
	c.val[1][0] ^= t[0];
	t = vextq_p64(t, c.val[1], 1);
	c.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) t);
	c.val[1] = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[1], 1);
	return (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) c.val[1]);
}

static inline poly64x2_t bf_red_lazy(poly64x2x2_t c) {
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

poly64x2_t bf_red_lazy1(poly64x2x2_t c);

poly64x2_t bf_red_from_lazy(poly64x2_t a);

static inline poly64x2_t bf_square(poly64x2_t a) {
	poly64x2_t c0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], a[0]));
	poly64x2_t c1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[1], a[1]));
	
	poly64x2_t t;
	t[0] = 0;
	t = vextq_p64(c1, t, 1);
	c1[0] ^= t[0];
	t = vextq_p64(t, c1, 1);
	c0 = (poly64x2_t) veorq_u64((uint64x2_t) c0, (uint64x2_t) t);
	c1 = (poly64x2_t) vshlq_n_u64((uint64x2_t) c1, 1);
	return (poly64x2_t) veorq_u64((uint64x2_t) c0, (uint64x2_t) c1);
}

static inline poly64x2_t bf_mull(poly64x2_t a, poly64x2_t b) {
	poly64x2_t t0, t1;
	poly64x2_t z = {0,0};

	poly64x2_t r0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], b[0]));
	poly64x2_t r1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a, b));
	t0 = vextq_p64(b, b, 1);
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], t0[0]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a, t0));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	r0 = (poly64x2_t) veorq_u64((uint64x2_t) r0, (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	r1 = (poly64x2_t) veorq_u64((uint64x2_t) r1, (uint64x2_t) t1);

	t0[0] = 0;
	t1[0] = r1[0] >> 63;
	t0 = vextq_p64(r1, t0, 1);
	r1[0] ^= t0[0];
	t1[0] ^= r1[0];
	t0 = vextq_p64(t0, t1, 1);
	r0 = (poly64x2_t) veorq_u64((uint64x2_t) r0, (uint64x2_t) t0);
	r1 = (poly64x2_t) vshlq_n_u64((uint64x2_t) r1, 1);
	return (poly64x2_t) veorq_u64((uint64x2_t) r0, (uint64x2_t) r1);
}

static inline poly64x2_t bf_multisquare_loop(poly64x2_t a, uint64_t n) {
	for(int i = 0; i < n; i++) {
		a = bf_red_psquare(bf_psquare(a));
	}
	return a;
}

poly64x2_t bf_inv(poly64x2_t a);

//Implementation alternatives:

poly64x2x2_t bf_pmull32(poly64x2_t a, poly64x2_t b);

poly64x2x2_t bf_pmull64(poly64x2_t a, poly64x2_t b);

poly64x2_t bf_red_generic(poly64x2x2_t c);

poly64x2_t bf_red_formula(poly64x2x2_t c);

poly64x2_t bf_red_neon(poly64x2x2_t c);

poly64x2_t bf_red_psquare_formula(poly64x2x2_t c);

poly64x2_t bf_red_psquare_neon(poly64x2x2_t c);

poly64x2_t bf_red_psquare_neonv2(poly64x2x2_t c);

poly64x2_t bf_red_lazy_formula(poly64x2x2_t c);

poly64x2_t bf_fermat_inv(poly64x2_t a);

poly64x2_t bf_addchain_inv(poly64x2_t a);

poly64x2_t bf_addchain_lookup_inv(poly64x2_t a);

void precomp_inv_table(uint64_t k);

void precomp_inv_tables();

poly64x2_t bf_multisquare_lookup_6(poly64x2_t a);

poly64x2_t bf_multisquare_lookup_12(poly64x2_t a);

poly64x2_t bf_multisquare_lookup_18(poly64x2_t a);

poly64x2_t bf_multisquare_lookup_30(poly64x2_t a);

poly64x2_t bf_multisquare_lookup_48(poly64x2_t a);

static inline poly64x2_t shift_right_carry(poly64x2_t a) {
	poly64x2_t c = a;
	c[0] /= 2;
	if(c[1] % 2 == 1)
		c[0] += pow2to63;
	c[1] /= 2;
	return c;
}

static inline uint64_t hasGreaterDeg(poly64x2_t u, poly64x2_t v) {
	uint64_t lz_u_h, lz_u_l, lz_v_h, lz_v_l;
	__asm (
				"CLZ %[lz_u_h], %[in_u_h]\n\t"
				"CLZ %[lz_u_l], %[in_u_l]\n\t"
				"CLZ %[lz_v_h], %[in_v_h]\n\t"
				"CLZ %[lz_v_l], %[in_v_l]\n\t"
				: [lz_u_h] "=r" (lz_u_h), [lz_u_l] "=r" (lz_u_l), [lz_v_h] "=r" (lz_v_h), [lz_v_l] "=r" (lz_v_l)  
				: [in_u_h] "r" (u[1]), [in_u_l] "r" (u[0]), [in_v_h] "r" (v[1]), [in_v_l] "r" (v[0])
			);
	uint64_t deg_u, deg_v;
	if(lz_u_h == 64) {
		deg_u = 63 - lz_u_l;
	} else {
		deg_u = 127 - lz_u_h;
	}
	if(lz_v_h == 64) {
		deg_v = 63 - lz_v_l;
	} else {
		deg_v = 127 - lz_v_h;
	}

	//printf("deg(u) = %lu deg(v) = %lu\n", deg_u, deg_v);
	return deg_u > deg_v;
}

#endif
