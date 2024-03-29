#include <stdio.h>

#include "ec.h"
#include "utils.h"
#include "ec_scalarmull.h"

ec_point_lproj ec_create_point_lproj(ef_intrl_elem x, ef_intrl_elem l, ef_intrl_elem z) {
	ec_point_lproj P;
	P.x = x;
	P.l = l;
	P.z = z;
	return P;
}

ec_point_laffine ec_create_point_laffine(ef_intrl_elem x, ef_intrl_elem l) {
	ec_point_laffine P;
	P.x = x;
	P.l = l;
	return P;
}

ec_point_lproj ec_laffine_to_lproj(ec_point_laffine P) {
	ec_point_lproj R;
	R.x = P.x;
	R.l = P.l;
	R.z = (ef_intrl_elem) {{{1, 0}, {0, 0}}};
	return R;
}

void ec_laffine_to_lproj_ptr(ec_point_laffine *P, ec_point_lproj *R) {
	R->x = P->x;
	R->l = P->l;
	R->z = (ef_intrl_elem) {{{1, 0}, {0, 0}}};
}

//Leads to undefined behavior for P == INFTY
ec_point_laffine ec_lproj_to_laffine(ec_point_lproj P, int in_const_time) {
	ef_intrl_elem Z_inv = ef_intrl_inv(P.z, in_const_time);
	ec_point_laffine R;
	R.x = ef_intrl_mull(P.x, Z_inv);
	R.l = ef_intrl_mull(P.l, Z_inv);
	return R;
}

void ec_lproj_to_laffine_ptr(ec_point_lproj *P, ec_point_laffine *R, int in_const_time) {
	ef_intrl_elem Z_inv = ef_intrl_inv(P->z, in_const_time);
	R->x = ef_intrl_mull(P->x, Z_inv);
	R->l = ef_intrl_mull(P->l, Z_inv);
}

void ec_print_expr(ec_point_lproj P) {
	printf("x: ");
	ef_intrl_print_expr_nl(P.x);
	printf(" l: ");
	ef_intrl_print_expr_nl(P.l);
	printf(" z: ");
	ef_intrl_print_expr_nl(P.z);
}

void ec_print_expr_laffine(ec_point_laffine P) {
	printf("x: ");
	ef_intrl_print_expr_nl(P.x);
	printf(" l: ");
	ef_intrl_print_expr_nl(P.l);
}

void ec_print_hex(ec_point_lproj P) {
	printf("x: ");
	ef_intrl_print_hex_nl(P.x);
	printf(" l: ");
	ef_intrl_print_hex_nl(P.l);
	printf(" z: ");
	ef_intrl_print_hex_nl(P.z);
}

void ec_print_hex_laffine(ec_point_laffine P) {
	printf("x: ");
	ef_intrl_print_hex_nl(P.x);
	printf(" l: ");
	ef_intrl_print_hex_nl(P.l);
}

uint64_t ec_is_on_curve(ec_point_lproj P) {
	ef_intrl_elem lhs = ef_intrl_mull(ef_intrl_add(ef_intrl_square(P.l), ef_intrl_add(ef_intrl_mull(P.l, P.z), ef_intrl_mull((ef_intrl_elem) A, ef_intrl_square(P.z)))), ef_intrl_square(P.x)); //(L^2 + LZ + AZ^2)X^2
	ef_intrl_elem rhs = ef_intrl_add(ef_intrl_square(ef_intrl_square(P.x)), ef_intrl_mull((ef_intrl_elem) B, ef_intrl_square(ef_intrl_square(P.z)))); //X^4 + BZ^4
	return ef_intrl_equal(lhs, rhs);
}

uint64_t ec_is_on_curve_laffine(ec_point_laffine P) {
	return ec_is_on_curve(ec_laffine_to_lproj(P));
}

uint64_t ec_equal_point_lproj(ec_point_lproj P, ec_point_lproj Q) {
	ef_intrl_elem zero = (ef_intrl_elem) {{{0,0}, {0,0}}};
	if(ec_is_on_curve(P) && ec_is_on_curve(Q) && ef_intrl_equal(P.z, zero) && ef_intrl_equal(Q.z, zero)) {
		return 1;
	}
	ec_point_laffine P_affine = ec_lproj_to_laffine(P, 0);
	ec_point_laffine Q_affine = ec_lproj_to_laffine(Q, 0);
	return ef_intrl_equal(P_affine.x, Q_affine.x) && ef_intrl_equal(P_affine.l, Q_affine.l);
}

uint64_t ec_equal_point_mixed(ec_point_laffine P, ec_point_lproj Q) {
	ef_intrl_elem zero = (ef_intrl_elem) {{{0,0}, {0,0}}};
	if(ef_intrl_equal(Q.z, zero)) {
		return 0;
	}
	ec_point_laffine Q_affine = ec_lproj_to_laffine(Q, 0);
	return ef_intrl_equal(P.x, Q_affine.x) && ef_intrl_equal(P.l, Q_affine.l);
}

uint64_t ec_equal_point_laffine(ec_point_laffine P, ec_point_laffine Q) {
	return ef_intrl_equal(P.x, Q.x) && ef_intrl_equal(P.l, Q.l);
}

// Generate random number in range [1, ORDER-1]
uint64x2x2_t ec_rand_scalar() {
	uint64x2x2_t order = (uint64x2x2_t) SUBGROUP_ORDER;
	uint64x2x2_t k;

	int in_range = 0;
	while (!in_range) {
		uint64_t a0 = rand_uint64();
		uint64_t a1 = rand_uint64();
		uint64_t a2 = rand_uint64();
		uint64_t a3 = rand_uint64();

		if (a3 > order.val[1][1]) continue;
		if (a3 == order.val[1][1] && a2 > order.val[1][0]) continue;
		if (a3 == order.val[1][1] && a2 == order.val[1][0] && a1 > order.val[0][1]) continue;
		if (a3 == order.val[1][1] && a2 == order.val[1][0] && a1 == order.val[0][1] && a0 >= order.val[0][0]) continue;

		in_range = 1;

		uint64x2_t p1 = { a0, a1 };
		uint64x2_t p2 = { a2, a3 };

		k.val[0] = p1;
		k.val[1] = p2;
	}

	return k;
}

ec_point_lproj ec_rand_point_lproj() {
	uint64x2x2_t k = ec_rand_scalar();

	return ec_scalarmull_single_lproj((ec_point_lproj) GEN, k);
}

ec_point_laffine ec_rand_point_laffine() {
	ec_point_lproj P = ec_rand_point_lproj();
	return ec_lproj_to_laffine(P, 0);
}

// Non constant implementation.
ec_point_lproj ec_add(ec_point_lproj P, ec_point_lproj Q) {
	if(ec_equal_point_lproj(P, (ec_point_lproj) INFTY)) {
		return Q;
	}
	if(ec_equal_point_lproj(Q, (ec_point_lproj) INFTY)) {
		return P;
	}
	if(ec_equal_point_lproj(P, ec_neg(Q))) {
		return (ec_point_lproj) INFTY;
	}
	if(ec_equal_point_lproj(P, Q)) {
		return ec_double(P);
	}

	return ec_add_unchecked(P, Q);
}

ec_point_lproj ec_add_unchecked(ec_point_lproj P, ec_point_lproj Q) {
	ef_intrl_elem u = ef_intrl_add(ef_intrl_mull(P.l, Q.z), ef_intrl_mull(Q.l, P.z)); // U = L_P * Z_Q + L_Q * Z_P
	ef_intrl_elem w1 = ef_intrl_mull(P.x, Q.z); //W1 = X_P * Z_Q
	ef_intrl_elem w2 = ef_intrl_mull(Q.x, P.z); //W2 = X_Q * Z_P
	ef_intrl_elem v = ef_intrl_square(ef_intrl_add(w1, w2)); //V = (X_P * Z_Q + X_Q * Z_P)^2
	ef_intrl_elem w3 = ef_intrl_mull(u, w2); //W3 = U * X_Q * Z_P
	ef_intrl_elem w4 = ef_intrl_mull(u, ef_intrl_mull(v, Q.z)); //W4 = U * V * Z_Q
	ec_point_lproj R;
	R.x = ef_intrl_mull(u, ef_intrl_mull(w1, w3));
	R.l = ef_intrl_add(ef_intrl_square(ef_intrl_add(w3, v)), ef_intrl_mull(w4, ef_intrl_add(P.l, P.z)));
	R.z = ef_intrl_mull(w4, P.z);
	return R;
}

void ec_add_ptr(ec_point_lproj *P, ec_point_lproj *Q, ec_point_lproj *R) {
	ef_intrl_elem zero = (ef_intrl_elem) {{{0,0}, {0,0}}};
	ec_point_lproj inf = (ec_point_lproj) INFTY;
	ec_point_lproj *tmp_res = &inf;

	uint64_t is_P_inf = ef_intrl_equal(P->z, zero);

	CSEL_PTR(is_P_inf, tmp_res, Q);
	
	uint64_t is_Q_inf = ef_intrl_equal(Q->z, zero);
	CSEL_PTR(is_Q_inf, tmp_res, P);

	ef_intrl_elem Xp = ef_intrl_mull(P->x, Q->z);
	ef_intrl_elem Lp = ef_intrl_mull(P->l, Q->z);
	ef_intrl_elem Xq = ef_intrl_mull(Q->x, P->z);
	ef_intrl_elem Lq = ef_intrl_mull(Q->l, P->z);
	ef_intrl_elem Z = ef_intrl_mull(P->z, Q->z);

	uint64_t equal_x = ef_intrl_equal(Xp, Xq);
	uint64_t equal_l = ef_intrl_equal(Lp, Lq);
	uint64_t neg_l = ef_intrl_equal(ef_intrl_add(Lp, Z), Lq);
	uint64_t neg = equal_x & neg_l;

	uint64_t equal = equal_x & equal_l;
	ec_point_lproj P2;
	ec_double_alt_ptr(P, &P2);
	CSEL_PTR(equal, tmp_res, &P2);
	
	uint64_t normal = (!is_P_inf) & (!is_Q_inf) & (!neg) & (!equal);
	ec_point_lproj PplusQ;
	ec_add_unchecked_ptr(P, Q, &PplusQ);
	CSEL_PTR(normal, tmp_res, &PplusQ);
	*R = *tmp_res;
}

void ec_add_unchecked_ptr(ec_point_lproj *P, ec_point_lproj *Q, ec_point_lproj *R) {
	ef_intrl_elem u = ef_intrl_add(ef_intrl_mull(P->l, Q->z), ef_intrl_mull(Q->l, P->z)); // U = L_P * Z_Q + L_Q * Z_P
	ef_intrl_elem w1 = ef_intrl_mull(P->x, Q->z); //W1 = X_P * Z_Q
	ef_intrl_elem w2 = ef_intrl_mull(Q->x, P->z); //W2 = X_Q * Z_P
	ef_intrl_elem v = ef_intrl_square(ef_intrl_add(w1, w2)); //V = (X_P * Z_Q + X_Q * Z_P)^2
	ef_intrl_elem w3 = ef_intrl_mull(u, w2); //W3 = U * X_Q * Z_P
	ef_intrl_elem w4 = ef_intrl_mull(u, ef_intrl_mull(v, Q->z)); //W4 = U * V * Z_Q
	R->x = ef_intrl_mull(u, ef_intrl_mull(w1, w3));
	R->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(w3, v)), ef_intrl_mull(w4, ef_intrl_add(P->l, P->z)));
	R->z = ef_intrl_mull(w4, P->z);
}

ec_point_lproj ec_add_mixed(ec_point_laffine P, ec_point_lproj Q) {
	if(ec_equal_point_lproj(Q, (ec_point_lproj) INFTY)) {
		return ec_laffine_to_lproj(P);
	}
	if(ec_equal_point_mixed(P, ec_neg(Q))) {
		return (ec_point_lproj) INFTY;
	}
	if(ec_equal_point_mixed(P, Q)) {
		return ec_double(Q);
	}

	return ec_add_mixed_unchecked(P, Q);
}

ec_point_lproj ec_add_mixed_unchecked(ec_point_laffine P, ec_point_lproj Q) {
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(P.l, Q.z), Q.l); //A = L_P * Z_Q + L_Q
	ef_intrl_elem F = ef_intrl_mull(P.x, Q.z); //X_P * Z_Q
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, Q.x)); //B = (X_P * Z_Q + X_Q)^2
	ef_intrl_elem H = ef_intrl_mull(E, Q.x); //A * Q.x
	//P.l.val[0] = bf_add(P.l.val[0], (poly64x2_t) {1, 0});
	ec_point_lproj R;
	R.x = ef_intrl_mull(ef_intrl_mull(E, F), H); //A * (X_P * Z_Q) * A * Q.x
	R.z = ef_intrl_mull(ef_intrl_mull(E, G), Q.z); //A * B * Z_Q
	R.l = ef_intrl_add(ef_intrl_square(ef_intrl_add(G, H)), ef_intrl_mull(R.z, ef_intrl_add(P.l, (ef_intrl_elem) {{{1, 0}, {0, 0}}}))); //(G+H)^2 + R.z + R.z * (P.l + 1)
	return R;
}

void ec_add_mixed_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *R) {
	ef_intrl_elem zero = (ef_intrl_elem) {{{0,0}, {0,0}}};
	ec_point_lproj inf = (ec_point_lproj) INFTY;
	ec_point_lproj *tmp_res = &inf;
	
	uint64_t is_Q_inf = ef_intrl_equal(Q->z, zero);
	ec_point_lproj Pproj;
	ec_laffine_to_lproj_ptr(P, &Pproj);
	CSEL_PTR(is_Q_inf, tmp_res, &Pproj);

	ef_intrl_elem Xp = ef_intrl_mull(P->x, Q->z);
	ef_intrl_elem Lp = ef_intrl_mull(P->l, Q->z);

	uint64_t equal_x = ef_intrl_equal(Xp, Q->x);
	uint64_t equal_l = ef_intrl_equal(Lp, Q->l);
	uint64_t neg_l = ef_intrl_equal(ef_intrl_add(Lp, Q->z), Q->l);
	uint64_t neg = equal_x & neg_l;

	uint64_t equal = equal_x & equal_l;
	ec_point_lproj P2;
	ec_double_mixed_ptr(P, &P2);
	CSEL_PTR(equal, tmp_res, &P2);
	
	uint64_t normal = (!is_Q_inf) & (!neg) & (!equal);
	ec_point_lproj PplusQ;
	ec_add_mixed_unchecked_ptr(P, Q, &PplusQ);
	CSEL_PTR(normal, tmp_res, &PplusQ);

	*R = *tmp_res;
}

void ec_add_mixed_unchecked_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *R) {
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(P->l, Q->z), Q->l); //A = L_P * Z_Q + L_Q
	ef_intrl_elem F = ef_intrl_mull(P->x, Q->z); //X_P * Z_Q
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, Q->x)); //B = (X_P * Z_Q + X_Q)^2
	ef_intrl_elem H = ef_intrl_mull(E, Q->x); //A * Q.x
	//P.l.val[0] = bf_add(P.l.val[0], (poly64x2_t) {1, 0});
	R->x = ef_intrl_mull(ef_intrl_mull(E, F), H); //A * (X_P * Z_Q) * A * Q.x
	R->z = ef_intrl_mull(ef_intrl_mull(E, G), Q->z); //A * B * Z_Q
	R->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(G, H)), ef_intrl_mull(R->z, ef_intrl_add(P->l, (ef_intrl_elem) {{{1, 0}, {0, 0}}}))); //(G+H)^2 + R.z + R.z * (P.l + 1)
}

ec_point_lproj ec_add_laffine_unchecked(ec_point_laffine P, ec_point_laffine Q) {
	ef_intrl_elem E = ef_intrl_add(P.l, Q.l); //A = L_P + L_Q
	ef_intrl_elem F = ef_intrl_square(ef_intrl_add(P.x, Q.x)); //B = (X_P + X_Q)^2
	ef_intrl_elem G = ef_intrl_mull(E, Q.x); //A * X_Q
	//P.l.val[0] = bf_add(P.l.val[0], (poly64x2_t) {1, 0});
	ec_point_lproj R;
	R.x = ef_intrl_mull(ef_intrl_mull(E, P.x), G); //A * X_P * A * X_Q
	R.z = ef_intrl_mull(E, F); //A * B * Z_Q
	R.l = ef_intrl_add(ef_intrl_square(ef_intrl_add(F, G)), ef_intrl_mull(R.z, ef_intrl_add(P.l, (ef_intrl_elem) {{{1, 0}, {0, 0}}}))); //(G+H)^2 + R.z + R.z * (P.l + 1)
	return R;
}

void ec_add_laffine_ptr(ec_point_laffine *P, ec_point_laffine *Q, ec_point_lproj *R) {
	ec_point_lproj inf = (ec_point_lproj) INFTY;
	ec_point_lproj *tmp_res = &inf;
	ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};

	uint64_t equal_x = ef_intrl_equal(P->x, Q->x);
	uint64_t equal_l = ef_intrl_equal(P->l, Q->l);
	uint64_t neg_l = ef_intrl_equal(ef_intrl_add(P->l, one), Q->l);
	uint64_t neg = equal_x & neg_l;

	uint64_t equal = equal_x & equal_l;
	ec_point_lproj P2;
	ec_double_mixed_ptr(P, &P2);
	CSEL_PTR(equal, tmp_res, &P2);
	
	uint64_t normal = (!neg) & (!equal);
	ec_point_lproj PplusQ;
	ec_add_laffine_unchecked_ptr(P, Q, &PplusQ);
	CSEL_PTR(normal, tmp_res, &PplusQ);

	*R = *tmp_res;
}

void ec_add_laffine_unchecked_ptr(ec_point_laffine *P, ec_point_laffine *Q, ec_point_lproj *R) {
	ef_intrl_elem E = ef_intrl_add(P->l, Q->l); //A = L_P + L_Q
	ef_intrl_elem F = ef_intrl_square(ef_intrl_add(P->x, Q->x)); //B = (X_P + X_Q)^2
	ef_intrl_elem G = ef_intrl_mull(E, Q->x); //A * X_Q
	//P.l.val[0] = bf_add(P.l.val[0], (poly64x2_t) {1, 0});
	R->x = ef_intrl_mull(ef_intrl_mull(E, P->x), G); //A * X_P * A * X_Q
	R->z = ef_intrl_mull(E, F); //A * B * Z_Q
	R->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(F, G)), ef_intrl_mull(R->z, ef_intrl_add(P->l, (ef_intrl_elem) {{{1, 0}, {0, 0}}}))); //(G+H)^2 + R.z + R.z * (P.l + 1)
}

//P+-Q
void ec_add_sub_laffine_unchecked_ptr(ec_point_laffine *P, ec_point_laffine *Q, ec_point_lproj *Radd, ec_point_lproj *Rsub) {
	ef_intrl_elem E = ef_intrl_add(P->l, Q->l); //A = L_P + L_Q
	ef_intrl_elem F = ef_intrl_square(ef_intrl_add(P->x, Q->x)); //B = (X_P + X_Q)^2
	ef_intrl_elem G = ef_intrl_mull(P->x, Q->x); //X_P * X_Q
	ef_intrl_elem lPplusOne = P->l;
	lPplusOne.val[0][0] ^= 1;
	Radd->x = ef_intrl_mull(ef_intrl_square(E), G); //A² * (X_P * X_Q)
	Radd->z = ef_intrl_mull(E, F);
	Radd->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(ef_intrl_mull(E, Q->x), F)), ef_intrl_mull(Radd->z, lPplusOne));
	Rsub->x = ef_intrl_add(Radd->x, G);
	Rsub->l = ef_intrl_add(ef_intrl_add(Radd->l, ef_intrl_square(Q->x)), ef_intrl_mull(F, lPplusOne));
	Rsub->z = ef_intrl_add(Radd->z, F);
}

void ec_add_sub_mixed_unchecked_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *Radd, ec_point_lproj *Rsub) {
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(P->l, Q->z), Q->l); //A = L_P * Z_Q + L_Q
	ef_intrl_elem xPzQ = ef_intrl_mull(P->x, Q->z);
	ef_intrl_elem F = ef_intrl_square(ef_intrl_add(xPzQ, Q->x));
	ef_intrl_elem G = ef_intrl_mull(xPzQ, Q->x);
	ef_intrl_elem lPplusOne = P->l;
	lPplusOne.val[0][0] ^= 1;
	Radd->x = ef_intrl_mull(ef_intrl_square(E), G);
	Radd->z = ef_intrl_mull(ef_intrl_mull(E, F), Q->z);
	Radd->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(ef_intrl_mull(E, Q->x), F)), ef_intrl_mull(Radd->z, lPplusOne));
	ef_intrl_elem zQzQ = ef_intrl_square(Q->z);
	Rsub->x = ef_intrl_add(Radd->x, ef_intrl_mull(G, zQzQ));
	ef_intrl_elem H = ef_intrl_mull(zQzQ, F);
	Rsub->z = ef_intrl_add(Radd->z, H);
	Rsub->l = ef_intrl_add(ef_intrl_add(Radd->l, ef_intrl_square(ef_intrl_mull(Q->x, Q->z))), ef_intrl_mull(H, lPplusOne));
}

//Computes P + psi(P)
void ec_add_endo_laffine_unchecked_ptr(ec_point_laffine *P, ec_point_lproj *R) {
	//E = l + l + l1+u = l1+ u = {l1, 1}
	//F = (x + x + x1)² = x1² = {x1²,0}
	//G = E*x = (l1+u)(x0+x1u) = {(x0l1+x1) + (x1l1+ x0+x1)u
	//H = E * x1 = {x1*l1, x1}
	//R->x = G * (G+H) = E * x * (E*x + E*x1) = E * x * E * (x+x1)
	//R->z = E*F = (l1+u)*x1² = x1²l1+x1²u = {x1² * l1, x1²}
	//R->l = (G+H+F)² + R->z *(l+1)
	poly64x2x2_t x_nonintrl = ef_intrl_disentangle(P->x);
	poly64x2x2_t l_nonintrl = ef_intrl_disentangle(P->l);
	poly64x2x2_t G_nonintrl;
	G_nonintrl.val[0] = bf_add(bf_red_lazy(bf_pmull(x_nonintrl.val[0], l_nonintrl.val[1])), x_nonintrl.val[1]);
	poly64x2_t x1l1 = bf_red_lazy(bf_pmull(x_nonintrl.val[1], l_nonintrl.val[1]));
	G_nonintrl.val[1] = bf_add(bf_add(x1l1, x_nonintrl.val[0]), x_nonintrl.val[1]);
	ef_intrl_elem G = ef_intrl_interleave(G_nonintrl);
	poly64x2x2_t H_nonintrl;
	H_nonintrl.val[0] = x1l1;
	H_nonintrl.val[1] = x_nonintrl.val[1];
	ef_intrl_elem H = ef_intrl_interleave(H_nonintrl);
	
	R->x = ef_intrl_mull(G, ef_intrl_add(G, H));
	poly64x2x2_t Rz_nonintrl;
	Rz_nonintrl.val[1] = bf_red_psquare(bf_psquare(x_nonintrl.val[1]));
	Rz_nonintrl.val[0] = bf_red_lazy(bf_pmull(Rz_nonintrl.val[1], l_nonintrl.val[1]));
	
	//G = G+F in interleaved form
	G.val[0][0] ^= Rz_nonintrl.val[1][0];
	G.val[1][0] ^= Rz_nonintrl.val[1][1];
	
	R->z = ef_intrl_interleave(Rz_nonintrl);
	R->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(G,H)), ef_intrl_mull(R->z, ef_intrl_add(P->l, (ef_intrl_elem) {{{1, 0}, {0, 0}}})));
}

ec_point_lproj ec_double(ec_point_lproj P) {
	ef_intrl_elem Z_sqr = ef_intrl_square(P.z);
	ef_intrl_elem S = ef_intrl_mull(P.l, P.z); //U = L_P * Z_P
	ef_intrl_elem T = ef_intrl_add(ef_intrl_square(P.l), ef_intrl_add(S, ef_intrl_mull_A(Z_sqr))); //T = L_P^2 + (L_P * Z_P) + A * Z_P^2
	ec_point_lproj R;
	R.x = ef_intrl_square(T);
	R.z = ef_intrl_mull(T, Z_sqr);
	R.l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_square(ef_intrl_mull(P.x, P.z)), R.x), ef_intrl_mull(T, S)), R.z);
	return R;
}

void ec_double_ptr(ec_point_lproj *P, ec_point_lproj *R) {
	ef_intrl_elem Z_sqr = ef_intrl_square(P->z);
	ef_intrl_elem S = ef_intrl_mull(P->l, P->z); //U = L_P * Z_P
	ef_intrl_elem T = ef_intrl_add(ef_intrl_square(P->l), ef_intrl_add(S, ef_intrl_mull_A(Z_sqr))); //T = L_P^2 + (L_P * Z_P) + A * Z_P^2
	R->x = ef_intrl_square(T);
	R->z = ef_intrl_mull(T, Z_sqr);
	R->l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_square(ef_intrl_mull(P->x, P->z)), R->x), ef_intrl_mull(T, S)), R->z);
}

ec_point_lproj ec_double_mixed(ec_point_laffine P) {
	ec_point_lproj R;
	R.z = ef_intrl_add(ef_intrl_square(P.l), ef_intrl_add(P.l, (ef_intrl_elem) A)); //T = L_P^2 + L_P + A
	R.x = ef_intrl_square(R.z);
	R.l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_square(P.x), R.x), ef_intrl_mull(R.z, P.l)), R.z);
	return R;
}

void ec_double_mixed_ptr(ec_point_laffine *P, ec_point_lproj *R) {
	R->z = ef_intrl_add(ef_intrl_square(P->l), ef_intrl_add(P->l, (ef_intrl_elem) A)); //T = L_P^2 + L_P + A
	R->x = ef_intrl_square(R->z);
	R->l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_square(P->x), R->x), ef_intrl_mull(R->z, P->l)), R->z);
}

void ec_double_alt_ptr(ec_point_lproj *P, ec_point_lproj *R) {
	ef_intrl_elem Z_sqr = ef_intrl_square(P->z);
	ef_intrl_elem Z_pow4 = ef_intrl_square(Z_sqr);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_square(P->l), ef_intrl_add(ef_intrl_mull(P->l, P->z), ef_intrl_mull_A(Z_sqr))); //W = L_P^2 + (L_P * Z_P) + A * Z_P^2
	ef_intrl_elem U = ef_intrl_square(ef_intrl_add(P->l, P->x));
	R->x = ef_intrl_square(T);
	R->z = ef_intrl_mull(T, Z_sqr);
	R->l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_mull(U, ef_intrl_add(ef_intrl_add(U, T), Z_sqr)), ef_intrl_mull_A(ef_intrl_mull_A(Z_pow4))), ef_intrl_mull_B(Z_pow4)), R->x), ef_intrl_mull_A(R->z)), R->z);
}

ec_point_lproj ec_double_then_add(ec_point_laffine P, ec_point_lproj Q) {
	ef_intrl_elem LQ_sqr = ef_intrl_square(Q.l);
	ef_intrl_elem ZQ_sqr = ef_intrl_square(Q.z);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_add(LQ_sqr, ef_intrl_mull(Q.l, Q.z)), ef_intrl_mull_A(ZQ_sqr));

	//ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};

	ef_intrl_elem tmp = P.l;
	tmp.val[0][0] ^= 1;
	tmp.val[0][1] ^= 1;
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(ef_intrl_square(Q.x), ZQ_sqr), ef_intrl_mull(T,ef_intrl_add(LQ_sqr, ef_intrl_mull(tmp, ZQ_sqr))));
	ef_intrl_elem F = ef_intrl_mull(P.x, ZQ_sqr);
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, T));
	tmp.val[0][1] ^= 1;
	ec_point_lproj R;
	R.x = ef_intrl_mull(F, ef_intrl_square(E));
	R.z = ef_intrl_mull(ef_intrl_mull(E, G), ZQ_sqr);
	R.l = ef_intrl_add(ef_intrl_mull(T, ef_intrl_square(ef_intrl_add(E, G))), ef_intrl_mull(tmp, R.z));
	return R;
}

void ec_double_then_add_ptr(ec_point_laffine *P, ec_point_lproj *Q, ec_point_lproj *R) {
	ef_intrl_elem LQ_sqr = ef_intrl_square(Q->l);
	ef_intrl_elem ZQ_sqr = ef_intrl_square(Q->z);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_add(LQ_sqr, ef_intrl_mull(Q->l, Q->z)), ef_intrl_mull_A(ZQ_sqr));

	//ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};

	ef_intrl_elem tmp = P->l;
	tmp.val[0][0] ^= 1;
	tmp.val[0][1] ^= 1;
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(ef_intrl_square(Q->x), ZQ_sqr), ef_intrl_mull(T,ef_intrl_add(LQ_sqr, ef_intrl_mull(tmp, ZQ_sqr))));
	ef_intrl_elem F = ef_intrl_mull(P->x, ZQ_sqr);
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, T));
	tmp.val[0][1] ^= 1;
	R->x = ef_intrl_mull(F, ef_intrl_square(E));
	R->z = ef_intrl_mull(ef_intrl_mull(E, G), ZQ_sqr);
	R->l = ef_intrl_add(ef_intrl_mull(T, ef_intrl_square(ef_intrl_add(E, G))), ef_intrl_mull(tmp, R->z));
}

ec_point_lproj ec_double_then_add_nonatomic(ec_point_laffine P, ec_point_lproj Q) {
	return ec_add_mixed_unchecked(P, ec_double(Q));
}

ec_point_lproj ec_double_then_addtwo(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q) {
	ef_intrl_elem LP1_plus_1 = P1.l;
	LP1_plus_1.val[0][0] ^= 1;
	ef_intrl_elem LQ_sqr = ef_intrl_square(Q.l);
	ef_intrl_elem ZQ_sqr = ef_intrl_square(Q.z);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_add(LQ_sqr, ef_intrl_mull(Q.l, Q.z)), ef_intrl_mull_A(ZQ_sqr)); //L_Q^2 + L_Q*Z_Q + A*Z_Q^2
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(ef_intrl_square(Q.x), ZQ_sqr), ef_intrl_mull(T, ef_intrl_add(LQ_sqr, ef_intrl_mull(ef_intrl_add((ef_intrl_elem) A, LP1_plus_1), ZQ_sqr)))); //X_Q^2*Z_Q^2 + T * (L_Q^2 + (a + L_P1 + 1)*Z_Q^2)
	ef_intrl_elem F = ef_intrl_mull(P1.x, ZQ_sqr); //X_P1 * Z_Q^2
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, T)); //(F+T)^2
	ef_intrl_elem H = ef_intrl_mull(ef_intrl_square(E), F); //X_2Q+P1 = F * E^2
	ef_intrl_elem I = ef_intrl_mull(ef_intrl_mull(E, G), ZQ_sqr); //Z_2Q+P1 = E * G * Z_Q^2
	ef_intrl_elem J = ef_intrl_add(ef_intrl_mull(ef_intrl_add(LP1_plus_1, P2.l), I), ef_intrl_mull(T, ef_intrl_square(ef_intrl_add(E, G)))); //(L_P1 + L_P2 + 1)*I+T*(E+G)^2
	ef_intrl_elem K = ef_intrl_mull(P2.x, I); // X_P2 * I
	ef_intrl_elem L = ef_intrl_square(ef_intrl_add(H, K)); //(H+K)^2
	ef_intrl_elem M = ef_intrl_mull(H, J); //H * J
	ec_point_lproj R;
	R.x = ef_intrl_mull(ef_intrl_mull(J, K), M);
	R.z = ef_intrl_mull(ef_intrl_mull(I, J), L);
	P2.l.val[0][0] ^= 1;
	R.l = ef_intrl_add(ef_intrl_square(ef_intrl_add(L, M)), ef_intrl_mull(R.z, P2.l));
	return R;
}

void ec_double_then_addtwo_ptr(ec_point_laffine *P1, ec_point_laffine *P2, ec_point_lproj *Q, ec_point_lproj *R) {
	ef_intrl_elem LP1_plus_1 = P1->l;
	LP1_plus_1.val[0][0] ^= 1;
	ef_intrl_elem LQ_sqr = ef_intrl_square(Q->l);
	ef_intrl_elem ZQ_sqr = ef_intrl_square(Q->z);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_add(LQ_sqr, ef_intrl_mull(Q->l, Q->z)), ef_intrl_mull_A(ZQ_sqr)); //L_Q^2 + L_Q*Z_Q + A*Z_Q^2
	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(ef_intrl_square(Q->x), ZQ_sqr), ef_intrl_mull(T, ef_intrl_add(LQ_sqr, ef_intrl_mull(ef_intrl_add((ef_intrl_elem) A, LP1_plus_1), ZQ_sqr)))); //X_Q^2*Z_Q^2 + T * (L_Q^2 + (a + L_P1 + 1)*Z_Q^2)
	ef_intrl_elem F = ef_intrl_mull(P1->x, ZQ_sqr); //X_P1 * Z_Q^2
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, T)); //(F+T)^2
	ef_intrl_elem H = ef_intrl_mull(ef_intrl_square(E), F); //X_2Q+P1 = F * E^2
	ef_intrl_elem I = ef_intrl_mull(ef_intrl_mull(E, G), ZQ_sqr); //Z_2Q+P1 = E * G * Z_Q^2
	ef_intrl_elem J = ef_intrl_add(ef_intrl_mull(ef_intrl_add(LP1_plus_1, P2->l), I), ef_intrl_mull(T, ef_intrl_square(ef_intrl_add(E, G)))); //(L_P1 + L_P2 + 1)*I+T*(E+G)^2
	ef_intrl_elem K = ef_intrl_mull(P2->x, I); // X_P2 * I
	ef_intrl_elem L = ef_intrl_square(ef_intrl_add(H, K)); //(H+K)^2
	ef_intrl_elem M = ef_intrl_mull(H, J); //H * J
	R->x = ef_intrl_mull(ef_intrl_mull(J, K), M);
	R->z = ef_intrl_mull(ef_intrl_mull(I, J), L);
	P2->l.val[0][0] ^= 1;
	R->l = ef_intrl_add(ef_intrl_square(ef_intrl_add(L, M)), ef_intrl_mull(R->z, P2->l));
}

ec_point_lproj ec_double_then_addtwo_nonatomic(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q) {
	return ec_add_mixed_unchecked(P2, ec_add_mixed_unchecked(P1, ec_double(Q)));
}

void ec_triple_mixed_ptr(ec_point_laffine *P, ec_point_lproj *R) {
	ef_intrl_elem T = ef_intrl_add(ef_intrl_square(P->l), ef_intrl_add(P->l, (ef_intrl_elem) A));
	ef_intrl_elem U = ef_intrl_add(ef_intrl_square(ef_intrl_add(P->x, T)), T);
	R->x = ef_intrl_mull(P->x, ef_intrl_square(U));
	R->z = ef_intrl_mull(ef_intrl_add(T, U), U);
	ef_intrl_elem lP_plus_1 = P->l;
	lP_plus_1.val[0][0] ^= 1;
	R->l = ef_intrl_add(ef_intrl_mull(ef_intrl_square(T), T), ef_intrl_mull(R->z, lP_plus_1));
}
