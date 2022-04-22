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

//Leads to undefined behavior for P == INFTY
ec_point_laffine ec_lproj_to_laffine(ec_point_lproj P) {
	ef_intrl_elem Z_inv = ef_intrl_inv(P.z);
	ec_point_laffine R;
	R.x = ef_intrl_mull(P.x, Z_inv);
	R.l = ef_intrl_mull(P.l, Z_inv);
	return R;
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
	ec_point_laffine P_affine = ec_lproj_to_laffine(P);
	ec_point_laffine Q_affine = ec_lproj_to_laffine(Q);
	return ef_intrl_equal(P_affine.x, Q_affine.x) && ef_intrl_equal(P_affine.l, Q_affine.l);
}

uint64_t ec_equal_point_mixed(ec_point_laffine P, ec_point_lproj Q) {
	ef_intrl_elem zero = (ef_intrl_elem) {{{0,0}, {0,0}}};
	if(ef_intrl_equal(Q.z, zero)) {
		return 0;
	}
	ec_point_laffine Q_affine = ec_lproj_to_laffine(Q);
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
	return ec_lproj_to_laffine(P);
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

/*
 * I think I have found a proof for why we don't need to check that P = -P here.
 * It seems that our projective lambda coords do not allow finite points of order 2:
 * 2*P = infty <-> P = -P <-> (x, l, z) = (x, l+z, z) => z=0
 *
 * For the case of the point at infty,
 * a quick calculation shows that this alg outputs infty again!
 */
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

ec_point_lproj ec_double_mixed(ec_point_laffine P) {
	ec_point_lproj R;
	R.z = ef_intrl_add(ef_intrl_square(P.l), ef_intrl_add(P.l, (ef_intrl_elem) A)); //T = L_P^2 + L_P + A
	R.x = ef_intrl_square(R.z);
	R.l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_intrl_square(P.x), R.x), ef_intrl_mull(R.z, P.l)), R.z);
	return R;
}

ec_point_lproj ec_double_alt(ec_point_lproj P) {
	ef_intrl_elem Z_sqr = ef_intrl_square(P.z);
	ef_intrl_elem Z_pow4 = ef_intrl_square(Z_sqr);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_square(P.l), ef_intrl_add(ef_intrl_mull(P.l, P.z), ef_intrl_mull_A(Z_sqr))); //W = L_P^2 + (L_P * Z_P) + A * Z_P^2
	ef_intrl_elem U = ef_intrl_square(ef_intrl_add(P.l, P.x));
	ec_point_lproj R;
	R.x = ef_intrl_square(T);
	R.z = ef_intrl_mull(T, Z_sqr);
	//R.l = ef_intrl_add(ef_intrl_add(ef_intrl_add(ef_mull(U, ef_add(ef_add(U, T), Z_sqr)), ef_mull_Apow2plusB(Z_pow4)), R.x), ef_mull_Aplus1(R.z));
	return R;
}

ec_point_lproj ec_double_then_add(ec_point_laffine P, ec_point_lproj Q) {
	ef_intrl_elem LQ_sqr = ef_intrl_square(Q.l);
	ef_intrl_elem ZQ_sqr = ef_intrl_square(Q.z);
	ef_intrl_elem T = ef_intrl_add(ef_intrl_add(LQ_sqr, ef_intrl_mull(Q.l, Q.z)), ef_intrl_mull_A(ZQ_sqr));

	ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};

	ef_intrl_elem E = ef_intrl_add(ef_intrl_mull(ef_intrl_square(Q.x), ZQ_sqr), ef_intrl_mull(T,ef_intrl_add(LQ_sqr, ef_intrl_mull(ef_intrl_add(ef_intrl_add((ef_intrl_elem)A, one), P.l), ZQ_sqr))));
	ef_intrl_elem F = ef_intrl_mull(P.x, ZQ_sqr);
	ef_intrl_elem G = ef_intrl_square(ef_intrl_add(F, T));

	ec_point_lproj R;
	R.x = ef_intrl_mull(F, ef_intrl_square(E));
	R.z = ef_intrl_mull(ef_intrl_mull(E, G), ZQ_sqr);
	R.l = ef_intrl_add(ef_intrl_mull(T, ef_intrl_square(ef_intrl_add(E, G))), ef_intrl_mull(ef_intrl_add(P.l, one), R.z));
	return R;
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

ec_point_lproj ec_double_then_addtwo_nonatomic(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q) {
	return ec_add_mixed_unchecked(P2, ec_add_mixed_unchecked(P1, ec_double(Q)));
}
