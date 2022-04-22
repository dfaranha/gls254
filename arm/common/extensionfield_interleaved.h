#include <arm_neon.h>
#include "extensionfield.h"

#ifndef EXTENSIONFIELD_INTERLEAVED_H
#define EXTENSIONFIELD_INTERLEAVED_H

typedef poly64x2x2_t ef_intrl_elem;

typedef struct {
	poly64x2_t val[4];
} ef_intrl_elem_unred;

void ef_intrl_print_expr(ef_intrl_elem a);

void ef_intrl_print_expr_nl(ef_intrl_elem a);

void ef_intrl_print_hex(ef_intrl_elem a);

void ef_intrl_print_hex_nl(ef_intrl_elem a);

void ef_intrl_print_unred_expr(ef_intrl_elem_unred a);

void ef_intrl_print_unred_expr_nl(ef_intrl_elem_unred a);

ef_elem ef_intrl_disentangle(ef_intrl_elem a);

ef_intrl_elem ef_intrl_interleave(ef_elem a);

poly64x2x2_t ef_intrl_disentangle_unred_lower(ef_intrl_elem_unred c);

poly64x2x2_t ef_intrl_disentangle_unred_higher(ef_intrl_elem_unred c);

ef_intrl_elem ef_intrl_rand_elem();

ef_intrl_elem_unred ef_intrl_rand_unred_elem();

uint64_t ef_intrl_equal(ef_intrl_elem a, ef_intrl_elem b);

static inline ef_intrl_elem ef_intrl_add(ef_intrl_elem a, ef_intrl_elem b) {
	ef_intrl_elem res;
	res.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a.val[0], (uint64x2_t) b.val[0]);
	res.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a.val[1], (uint64x2_t) b.val[1]);
	return res;
}

//Reduces mod z*f(z)
static inline ef_intrl_elem ef_intrl_red(ef_intrl_elem_unred c) {
	poly64x2_t t;
	c.val[2] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[2], (uint64x2_t) c.val[3]);
	//t = (poly64x2_t) vshrq_n_u64((uint64x2_t) c.val[3], 63); //Will this ever be necessary?? (z^127)^2 = 2^254 so last bit will never be set.
	//c.val[2] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[2], (uint64x2_t) t);
	t = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[3], 1);
	c.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[1], (uint64x2_t) t);
	c.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[1], (uint64x2_t) c.val[2]);
	t = (poly64x2_t) vshrq_n_u64((uint64x2_t) c.val[2], 63);
	c.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[1], (uint64x2_t) t);
	t = (poly64x2_t) vshlq_n_u64((uint64x2_t) c.val[2], 1);
	c.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) t);
	return (ef_intrl_elem) {{c.val[0], c.val[1]}};
}

ef_intrl_elem ef_intrl_red_from_lazy(ef_intrl_elem c);

static inline ef_intrl_elem ef_intrl_square(ef_intrl_elem a) {
	ef_intrl_elem_unred c;
	poly64x2_t t0, t1, t2, t3;

	//a_1^2
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a.val[0], a.val[0]));
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a.val[1], a.val[1]));

	//a_0^2 + a_1^2
	t2 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], a.val[0][0]));
	t2 = (poly64x2_t) veorq_u64((uint64x2_t) t2, (uint64x2_t) t0);
	t3 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], a.val[1][0]));
	t3 = (poly64x2_t) veorq_u64((uint64x2_t) t3, (uint64x2_t) t1);

	//Combine
	c.val[0] = (poly64x2_t) vzip1q_u64((uint64x2_t) t2, (uint64x2_t) t0);
	c.val[1] = (poly64x2_t) vzip2q_u64((uint64x2_t) t2, (uint64x2_t) t0);
	c.val[2] = (poly64x2_t) vzip1q_u64((uint64x2_t) t3, (uint64x2_t) t1);
	c.val[3] = (poly64x2_t) vzip2q_u64((uint64x2_t) t3, (uint64x2_t) t1);

	return ef_intrl_red(c);
}

static inline ef_intrl_elem ef_intrl_mull_A(ef_intrl_elem a) {
	poly64x2_t t = vextq_p64(a.val[0], a.val[0], 1);
	a.val[0][0] = 0;
	a.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a.val[0], (uint64x2_t) t);
	t = vextq_p64(a.val[1], a.val[1], 1);
	a.val[1][0] = 0;
	a.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a.val[1], (uint64x2_t) t);
	return a;
}

static inline ef_intrl_elem ef_intrl_mull(ef_intrl_elem a, ef_intrl_elem b) {
	ef_intrl_elem_unred c;
	poly64x2_t t0, t1;
	poly64x2_t z = vdupq_n_p64(0);
	poly64x2x2_t a0b0, a1b1, bigprod;

	//a0 * b0
	a0b0.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], b.val[0][0]));
	a0b0.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], b.val[1][0]));
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], b.val[1][0]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], b.val[0][0]));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	a0b0.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a0b0.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	a0b0.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a0b0.val[1], (uint64x2_t) t1);

	//a1 * b1
	a1b1.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a.val[0], b.val[0]));
	a1b1.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a.val[1], b.val[1]));
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a.val[0], b.val[1]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(a.val[1], b.val[0]));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	a1b1.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a1b1.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	a1b1.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a1b1.val[1], (uint64x2_t) t1);

	//a0+a1
	t0 = vextq_p64(a.val[0], z, 1);
	a.val[0][0] ^= t0[0];
	t0 = vextq_p64(a.val[1], z, 1);
	a.val[1][0] ^= t0[0];

	//b0+b1
	t0 = vextq_p64(b.val[0], z, 1);
	b.val[0][0] ^= t0[0];
	t0 = vextq_p64(b.val[1], z, 1);
	b.val[1][0] ^= t0[0];

	//(a0+a1)*(b0+b1)
	bigprod.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], b.val[0][0]));
	bigprod.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], b.val[1][0]));
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], b.val[1][0]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], b.val[0][0]));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	bigprod.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) bigprod.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	bigprod.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) bigprod.val[1], (uint64x2_t) t1);

	//Combine
	a1b1.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a0b0.val[0], (uint64x2_t) a1b1.val[0]);
	a1b1.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a0b0.val[1], (uint64x2_t) a1b1.val[1]);
	bigprod.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a0b0.val[0], (uint64x2_t) bigprod.val[0]);
	bigprod.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a0b0.val[1], (uint64x2_t) bigprod.val[1]);

	c.val[0] = (poly64x2_t) vzip1q_u64((uint64x2_t) a1b1.val[0], (uint64x2_t) bigprod.val[0]);
	c.val[1] = (poly64x2_t) vzip2q_u64((uint64x2_t) a1b1.val[0], (uint64x2_t) bigprod.val[0]);
	c.val[2] = (poly64x2_t) vzip1q_u64((uint64x2_t) a1b1.val[1], (uint64x2_t) bigprod.val[1]);
	c.val[3] = (poly64x2_t) vzip2q_u64((uint64x2_t) a1b1.val[1], (uint64x2_t) bigprod.val[1]);

	return ef_intrl_red(c);
}

ef_intrl_elem ef_intrl_inv(ef_intrl_elem a);

void ef_intrl_sim_inv(ef_intrl_elem inputs[], ef_intrl_elem outputs[], uint64_t len);

#endif
