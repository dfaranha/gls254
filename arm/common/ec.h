#include "basefield.h"
#include "extensionfield_interleaved.h"

#ifndef EC_H
#define EC_H

#define A {{{0, 1}, {0, 0}}}
#define B {{{134217729, 0}, {0,0}}} //z^27 + 1

#define GEN {{{{0x5BE5F4EB93D8712A, 0x9D1932CB5FA5B9BF}, {0x47E70D2DCA8C7210, 0x25F2F29FCBDEC78E}}}, {{{0x97FBBBBFEB3A8AB4, 0x25BE90C01E0E9B06}, {0x1A1764D658204447, 0x0B3834B048C217C1}}} , {{{1,0}, {0,0}}}}
#define INFTY {{{{1, 0}, {0, 0}}}, {{{1, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}}}

#define SUBGROUP_ORDER {{{0x877DABA2A44750A5, 0xDAC40D1195270779},{0xFFFFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFFFF}}}

typedef struct {
	ef_intrl_elem x, l, z;
} ec_point_lproj;

typedef struct {
	ef_intrl_elem x, l;
} ec_point_laffine;

ec_point_lproj ec_create_point_lproj(ef_intrl_elem x, ef_intrl_elem l, ef_intrl_elem z);

ec_point_laffine ec_create_point_laffine(ef_intrl_elem x, ef_intrl_elem l);

void ec_print_expr(ec_point_lproj P);

void ec_print_expr_laffine(ec_point_laffine P);

void ec_print_hex(ec_point_lproj P);

void ec_print_hex_laffine(ec_point_laffine P);

ec_point_lproj ec_laffine_to_lproj(ec_point_laffine P);

void ec_laffine_to_lproj_ptr(ec_point_laffine *P, ec_point_lproj *R);

ec_point_laffine ec_lproj_to_laffine(ec_point_lproj P, int in_const_time);

void ec_lproj_to_laffine_ptr(ec_point_lproj *P, ec_point_laffine *R, int in_const_time);

uint64_t ec_is_on_curve(ec_point_lproj P);

uint64_t ec_is_on_curve_laffine(ec_point_laffine P);

uint64_t ec_equal_point_lproj(ec_point_lproj P, ec_point_lproj Q);

uint64_t ec_equal_point_mixed(ec_point_laffine P, ec_point_lproj Q);

uint64_t ec_equal_point_laffine(ec_point_laffine P, ec_point_laffine Q);

uint64x2x2_t ec_rand_scalar();

ec_point_lproj ec_rand_point_lproj();

ec_point_laffine ec_rand_point_laffine();

static inline ec_point_lproj ec_neg(ec_point_lproj P) {
	P.l = ef_intrl_add(P.l, P.z);
	return P;
}

static inline void ec_neg_mut(ec_point_lproj *P) {
	P->l = ef_intrl_add(P->l, P->z);
}

static inline ec_point_laffine ec_neg_laffine(ec_point_laffine P) {
	P.l.val[0][0] ^= 1;
	return P;
}

static inline void ec_neg_laffine_mut(ec_point_laffine *P) {
	P->l.val[0][0] ^= 1;
}

ec_point_lproj ec_add(ec_point_lproj P1, ec_point_lproj P2);

ec_point_lproj ec_add_unchecked(ec_point_lproj P1, ec_point_lproj P2);

void ec_add_ptr(ec_point_lproj *P, ec_point_lproj *Q, ec_point_lproj *R);

void ec_add_unchecked_ptr(ec_point_lproj *P, ec_point_lproj *Q, ec_point_lproj *R);

ec_point_lproj ec_add_mixed(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_add_mixed_unchecked(ec_point_laffine P, ec_point_lproj Q);

void ec_add_mixed_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *R);

void ec_add_mixed_unchecked_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *R);

ec_point_lproj ec_add_laffine_unchecked(ec_point_laffine P, ec_point_laffine Q);

void ec_add_laffine_unchecked_ptr(ec_point_laffine *P, ec_point_laffine *Q, ec_point_lproj *R);

void ec_add_laffine_ptr(ec_point_laffine *P, ec_point_laffine *Q, ec_point_lproj *R);

void ec_add_sub_laffine_unchecked_ptr(ec_point_laffine *P, ec_point_laffine *Q, ec_point_lproj *Radd, ec_point_lproj *Rsub);

void ec_add_sub_mixed_unchecked_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *Radd, ec_point_lproj *Rsub);

void ec_add_endo_laffine_unchecked_ptr(ec_point_laffine *P, ec_point_lproj *R);

ec_point_lproj ec_double(ec_point_lproj P);

void ec_double_ptr(ec_point_lproj *P, ec_point_lproj *R);

ec_point_lproj ec_double_mixed(ec_point_laffine P);

void ec_double_mixed_ptr(ec_point_laffine *P, ec_point_lproj *R);

void ec_double_alt_ptr(ec_point_lproj *P, ec_point_lproj *R);

ec_point_lproj ec_double_then_add(ec_point_laffine P, ec_point_lproj Q);

void ec_double_then_add_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *R);

ec_point_lproj ec_double_then_add_nonatomic(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_double_then_addtwo(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q);

void ec_double_then_addtwo_ptr(ec_point_laffine *P1, ec_point_laffine *P2, ec_point_lproj *Q, ec_point_lproj *R);

ec_point_lproj ec_double_then_addtwo_nonatomic(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q);

void ec_triple_mixed_ptr(ec_point_laffine *P, ec_point_lproj *R);

static inline ec_point_laffine ec_endo_laffine(ec_point_laffine P) {
	/*P.x.val[0] = bf_add(P.x.val[0], P.x.val[1]);
	P.l.val[0] = bf_add(P.l.val[0], P.l.val[1]);
	P.l.val[1] = bf_add(P.l.val[1], (poly64x2_t) {1,0});*/

	//x_1u + (x0+x1)
	poly64x2_t t = vextq_p64(P.x.val[0], P.x.val[0], 1);
	P.x.val[0][0] ^= t[0];
	t = vextq_p64(P.x.val[1], P.x.val[1], 1);
	P.x.val[1][0] ^= t[0];

	//(l_1+1)u + (l0+l1)
	t[0] = 1;
	t = vextq_p64(P.l.val[0], t, 1);
	P.l.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) P.l.val[0], (uint64x2_t) t);
	t = vextq_p64(P.l.val[1], P.l.val[1], 1);
	P.l.val[1][0] ^= t[0];

	return P;
}

static inline void ec_endo_laffine_ptr(ec_point_laffine *P, ec_point_laffine *R) {
	/*P.x.val[0] = bf_add(P.x.val[0], P.x.val[1]);
	P.l.val[0] = bf_add(P.l.val[0], P.l.val[1]);
	P.l.val[1] = bf_add(P.l.val[1], (poly64x2_t) {1,0});*/

	//x_1u + (x0+x1)
	poly64x2_t t = vextq_p64(P->x.val[0], P->x.val[0], 1);
	R->x = P->x;
	R->x.val[0][0] ^= t[0];
	t = vextq_p64(P->x.val[1], P->x.val[1], 1);
	R->x.val[1][0] ^= t[0];

	//(l_1+1)u + (l0+l1)
	R->l = P->l;
	t[0] = 1;
	t = vextq_p64(P->l.val[0], t, 1);
	R->l.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) P->l.val[0], (uint64x2_t) t);
	t = vextq_p64(P->l.val[1], P->l.val[1], 1);
	R->l.val[1][0] ^= t[0];
}

static inline ec_point_lproj ec_endo_lproj(ec_point_lproj P) {
	ec_point_lproj R;
	poly64x2_t t = vextq_p64(P.x.val[0], P.x.val[0], 1);
	R.x = P.x;
	R.x.val[0][0] ^= t[0];
	t = vextq_p64(P.x.val[1], P.x.val[1], 1);
	R.x.val[1][0] ^= t[0];

	R.l.val[0] = vextq_p64(P.z.val[0], P.z.val[0], 1);
	R.z = P.z;
	R.z.val[0][0] ^= R.l.val[0][0];
	R.l.val[1] = vextq_p64(P.z.val[1], P.z.val[1], 1);
	R.z.val[1][0] ^= R.l.val[1][0];

	t = vextq_p64(P.l.val[0], P.l.val[0], 1);
	R.l = ef_intrl_add(R.l, P.l);
	R.l.val[0][0] ^= t[0];
	t = vextq_p64(P.l.val[1], P.l.val[1], 1);
	R.l.val[1][0] ^= t[0];
	
	return R;
}

#endif
