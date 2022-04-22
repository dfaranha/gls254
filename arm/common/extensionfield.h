#include <arm_neon.h>
#include "basefield.h"

#ifndef EXTENSIONFIELD_H
#define EXTENSIONFIELD_H

typedef poly64x2x2_t ef_elem;

static inline ef_elem ef_create_elem(poly64x2_t p0, poly64x2_t p1) {
	ef_elem a;
	a.val[0] = p0;
	a.val[1] = p1;
	return a;
}

void ef_print_expr(ef_elem a);

void ef_print_expr_nl(ef_elem a);

void ef_print_hex(ef_elem a);

void ef_print_hex_nl(ef_elem a);

ef_elem ef_rand_elem();

uint64_t ef_equal(ef_elem a, ef_elem b);

static inline ef_elem ef_add(ef_elem a, ef_elem b) {
	return ef_create_elem(bf_add(a.val[0], b.val[0]), bf_add(a.val[1], b.val[1]));
}

ef_elem ef_mull(ef_elem a, ef_elem b);

static inline ef_elem ef_mull_A(ef_elem a) {
	ef_elem r;
	r.val[0] = a.val[1];
	r.val[1] = bf_add(a.val[0], a.val[1]);
	return r;
}

static inline ef_elem ef_mull_Aplus1(ef_elem a) {
	ef_elem r;
	r.val[0] = bf_add(a.val[0], a.val[1]);
	r.val[1] = a.val[0];
	return r;
}

static inline ef_elem ef_mull_B(ef_elem a) {
	poly64x2_t t1, t2;
	ef_elem r;
	//r.val[0] *x^27+1
	t1 = (poly64x2_t) vshrq_n_u64((uint64x2_t) a.val[0], 37);
	t2 = vextq_p64(t1, t1, 1);
	t2[0] <<= 1;
	t1[0] = 0;
	t1 = (poly64x2_t) veorq_u64((uint64x2_t) t1, (uint64x2_t) t2);
	t2 = (poly64x2_t) vshlq_n_u64((uint64x2_t) a.val[0], 27);
	t1 = (poly64x2_t) veorq_u64((uint64x2_t) t1, (uint64x2_t) t2);
	r.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) a.val[0], (uint64x2_t) t1);
	//r.val[1] * (x^27 + 1)
	t1 = (poly64x2_t) vshrq_n_u64((uint64x2_t) a.val[1], 37);
	t2 = vextq_p64(t1, t1, 1);
	t2[0] <<= 1;
	t1[0] = 0;
	t1 = (poly64x2_t) veorq_u64((uint64x2_t) t1, (uint64x2_t) t2);
	t2 = (poly64x2_t) vshlq_n_u64((uint64x2_t) a.val[1], 27);
	t1 = (poly64x2_t) veorq_u64((uint64x2_t) t1, (uint64x2_t) t2);
	r.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) a.val[1], (uint64x2_t) t1);
	return r;
}

ef_elem ef_mull_Apow2plusB(ef_elem a);

ef_elem ef_square(ef_elem a);

ef_elem ef_inv(ef_elem a);

void ef_sim_inv(ef_elem inputs[], ef_elem outputs[], uint64_t len);

#endif
