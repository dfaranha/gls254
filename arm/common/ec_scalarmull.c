#include "ec_scalarmull.h"
#include <stdlib.h>
#include <string.h>

#define pow2to63 9223372036854775808U
#define pow2to64 {{{0, 1}, {0,0}}}

// Algorithm 3.27
ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k) {
	ec_point_lproj Q = (ec_point_lproj) INFTY;

	for(int i = 1; i >= 0; i--) {
		for(int j = 1; j >= 0; j--) {
			poly64_t c = pow2to63;

			for(int t = 63; t >= 0; t--) {
				Q = ec_double(Q);
				ec_point_lproj temp = ec_add_mixed(P, Q);

				uint64_t r0, val;
				asm ("TST %[k], %[c];"
					"CSEL %[r0], %[temp], %[q], NE;"
					: [r0] "=r" ( r0 )
					: [k] "r" (k.val[i][j]), [c] "r" (c), [q] "r" (&Q), [temp] "r" (&temp)
					);
				ec_point_lproj* ref = (ec_point_lproj*)r0;
				Q = *ref;

				c >>= 1;
			}
		}
	}
	return Q;
}

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k) {
  ec_point_lproj Q = (ec_point_lproj) INFTY;

  for(int i = 1; i >= 0; i--) {
    for(int j = 1; j >= 0; j--) {
      poly64_t c = pow2to63;

      for(int t = 63; t >= 0; t--) {
        Q = ec_double(Q);

        if (k.val[i][j] & c) {
          Q = ec_add(P, Q);
        }

        c >>= 1;
      }
    }
  }

  return Q;
}

// Algorithm 3.48
ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l) {
  ec_point_lproj R = (ec_point_lproj) INFTY;

  for(int i = 1; i >= 0; i--) {
    for(int j = 1; j >= 0; j--) {
      R = ec_scalarmull_single_lproj(R, (uint64x2x2_t) pow2to64);

      ec_point_lproj k_i_P = ec_scalarmull_single_lproj(P, (uint64x2x2_t) {{{k.val[i][j], 0}, {0,0}}});
      ec_point_lproj l_i_Q = ec_scalarmull_single_lproj(Q, (uint64x2x2_t) {{{l.val[i][j], 0}, {0,0}}});

      R = ec_add(R, ec_add(k_i_P, l_i_Q));
    }
  }

  return R;
}

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k) {
	ec_point_lproj new;
	uint64x1_t old_ptr, new_ptr, tmp, digitval;

	ec_split_scalar decomp = ec_scalar_decomp(k);
	ec_point_laffine Q = ec_endo_laffine(P);

	digitval[0]=1;
	ec_point_laffine P_neg = ec_neg_laffine(P);
	CMOV(tmp, decomp.k1_sign, digitval, P, P_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine Q_neg = ec_neg_laffine(Q);
	CMOV(tmp, decomp.k2_sign, digitval, Q, Q_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj R = (ec_point_lproj) INFTY;

	for(int i = 1; i >= 0; i--) {
		digitval[0] = pow2to63;
		for(int digit = 63; digit >= 0; digit--) {
			R = ec_double(R);

			new = ec_add_mixed(P, R);
			CMOV(tmp, decomp.k1[i], digitval, R, new, old_ptr, new_ptr, typeof(ec_point_lproj));

			new = ec_add_mixed(Q, R);
			CMOV(tmp, decomp.k2[i], digitval, R, new, old_ptr, new_ptr, typeof(ec_point_lproj));

			digitval[0] >>= 1;
		}
	}
	return ec_lproj_to_laffine(R, 1);
}

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 64;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 3, rec_k1, l);
	reg_rec(decomp.k2, 3, rec_k2, l);

	// Precomputation
	ec_point_laffine table[2];
	precompute_w3(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(-k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_mixed(P1, ec_laffine_to_lproj(P2));

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(Q);

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(-k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	return ec_lproj_to_laffine(Q, 1);
}

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 43;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l);
	reg_rec(decomp.k2, 4, rec_k2, l);

	// Precomputation
	ec_point_laffine table[4];
	precompute_w4(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(-k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	lin_pass_w4(&P1, &P2, &table, k1_val/2, k2_val/2);
	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(Q));

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(-k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w4(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	return ec_lproj_to_laffine(Q, 1);
}

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 32;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	// Compute recodings
	signed char rec_k1[l], rec_k2[l];

	reg_rec(decomp.k1, 5, rec_k1, l);
	reg_rec(decomp.k2, 5, rec_k2, l);

	// Precomputation
	ec_point_laffine table[8];
	precompute_w5(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(-k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1, P2;
	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(ec_double(Q)));

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(-k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		//P1 = table[k1_val/2];
		//P2 = table[k2_val/2];
		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	return ec_lproj_to_laffine(Q, 1);
}

void ec_scalarmull_single_endo_w5_randaccess_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R) {
	uint64_t new_ptr, con = 1;
	int l = 32;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	// Compute recodings
	signed char rec_k1[l], rec_k2[l];

	reg_rec(decomp.k1, 5, rec_k1, l);
	reg_rec(decomp.k2, 5, rec_k2, l);

	// Precomputation
	ec_point_laffine table[8];
	precompute_w5_ptr(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(- k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(-k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1, P2;
	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q, Q1, Q2;
	ec_add_laffine_unchecked_ptr(&P1, &P2, &Q);

	for(int i=l-2; i>=0; i--) {
		ec_double_alt_ptr(&Q, &Q1);
		ec_double_alt_ptr(&Q1, &Q2);
		ec_double_alt_ptr(&Q2, &Q1);

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(- k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(- k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));
		if(i == 0) {
			ec_double_alt_ptr(&Q1, &Q2);
			ec_add_mixed_ptr(&P1, &Q2, &Q1);
			ec_add_mixed_ptr(&P2, &Q1, &Q);
		} else {
			ec_double_then_addtwo_ptr(&P1, &P2, &Q1, &Q);
		}
	}

	ec_lproj_to_laffine_ptr(&Q, R, 1);
}

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 26;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 6, rec_k1, l);
	reg_rec(decomp.k2, 6, rec_k2, l);

	// Precomputation
	ec_point_laffine table[16];
	precompute_w6(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(- k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(- k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	lin_pass_w6(&P1, &P2, table, k1_val/2, k2_val/2);
	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));


	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);
	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(ec_double(ec_double(Q))));

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(- k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(- k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w6(&P1, &P2, table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		//Negate p1 by k1_sign
		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
		//Q = ec_double_then_add(P1, Q);
		//Q = ec_add_mixed_unchecked(P2, Q);
	}

	return ec_lproj_to_laffine(Q, 1);
}

void precompute_w3(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj Pl = ec_laffine_to_lproj(P);
	ec_point_lproj P3 = ec_double_then_add(P, Pl);

	ef_intrl_elem P3Z_inv = ef_intrl_inv(P3.z, 0);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, P3Z_inv), ef_intrl_mull(P3.l, P3Z_inv)};
}

void precompute_w4(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P5 = ec_double_then_add(P, P2);
	ec_point_lproj P7 = ec_double_then_add(P, P3);

	ef_intrl_elem inv_inputs[3] = {P3.z, P5.z, P7.z};
	ef_intrl_elem inv_outputs[3];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3, 0);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
}

void precompute_w5(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P4 = ec_double(P2);
	ec_point_lproj P5 = ec_add_mixed_unchecked(P, P4);
	ec_point_lproj P6 = ec_double(P3);
	ec_point_lproj P7 = ec_add_mixed_unchecked(P, P6);
	ec_point_lproj P9 = ec_double_then_add(P, P4);
	ec_point_lproj P11 = ec_double_then_add(P, P5);
	ec_point_lproj P13 = ec_double_then_add(P, P6);
	ec_point_lproj P15 = ec_double_then_add(P, P7);

	ef_intrl_elem inv_inputs[7] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z};
	ef_intrl_elem inv_outputs[7];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 7, 0);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
}

void precompute_w5_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ec_point_lproj P2;
	ec_double_mixed_ptr(P, &P2);
	ec_point_lproj P3;
	ec_add_mixed_unchecked_ptr(P, &P2, &P3);
	ec_point_lproj P4;
	ec_double_alt_ptr(&P2, &P4);
	ec_point_lproj P5;
	ec_add_mixed_unchecked_ptr(P, &P4, &P5);
	ec_point_lproj P6; 
	ec_double_alt_ptr(&P3, &P6);
	ec_point_lproj P7;
	ec_add_mixed_unchecked_ptr(P, &P6, &P7);
	ec_point_lproj P9;
	ec_double_then_add_ptr(P, &P4, &P9);
	ec_point_lproj P11;
	ec_double_then_add_ptr(P, &P5, &P11);
	ec_point_lproj P13;
	ec_double_then_add_ptr(P, &P6, &P13);
	ec_point_lproj P15;
	ec_double_then_add_ptr(P, &P7, &P15);

	ef_intrl_elem inv_inputs[7] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z};
	ef_intrl_elem inv_outputs[7];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 7, 0);

	table[0] = *P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
}

void precompute_w5_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ec_point_lproj tmp;
	ec_point_lproj table_proj[7];

	ec_triple_mixed_ptr(P, &table_proj[0]);
	for (int i = 0; i < 3; i++) {
		ec_double_alt_ptr(&table_proj[i], &tmp);
		ec_add_sub_mixed_unchecked_ptr(P, &tmp, &table_proj[2*i+2], &table_proj[2*i+1]);
		ec_neg_mut(&table_proj[2*i+1]);
	}
	
	ef_intrl_elem inv_inputs[7];
	ef_intrl_elem inv_outputs[7];
	for (int i = 0; i < 7; i++) {
		inv_inputs[i] = table_proj[i].z;
	}
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 7, 0);
	table[0] = *P;
	for (int i = 0; i < 7; i++) {
		table[i+1] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}
}

void precompute_w6(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P4 = ec_double(P2);
	ec_point_lproj P5 = ec_add_mixed_unchecked(P, P4);
	ec_point_lproj P6 = ec_double(P3);
	ec_point_lproj P7 = ec_add_mixed_unchecked(P, P6);
	ec_point_lproj P8 = ec_double(P4);
	ec_point_lproj P9 = ec_add_mixed_unchecked(P, P8);
	ec_point_lproj P10 = ec_double(P5);
	ec_point_lproj P11 = ec_add_mixed_unchecked(P, P10);
	ec_point_lproj P12 = ec_double(P6);
	ec_point_lproj P13 = ec_add_mixed_unchecked(P, P12);
	ec_point_lproj P14 = ec_double(P7);
	ec_point_lproj P15 = ec_add_mixed_unchecked(P, P14);
	ec_point_lproj P17 = ec_double_then_add(P, P8);
	ec_point_lproj P19 = ec_double_then_add(P, P9);
	ec_point_lproj P21 = ec_double_then_add(P, P10);
	ec_point_lproj P23 = ec_double_then_add(P, P11);
	ec_point_lproj P25 = ec_double_then_add(P, P12);
	ec_point_lproj P27 = ec_double_then_add(P, P13);
	ec_point_lproj P29 = ec_double_then_add(P, P14);
	ec_point_lproj P31 = ec_double_then_add(P, P15);

	ef_intrl_elem inv_inputs[15] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z, P17.z, P19.z, P21.z, P23.z, P25.z, P27.z, P29.z, P31.z};
	ef_intrl_elem inv_outputs[15];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 15, 0);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
	table[8] = (ec_point_laffine) {ef_intrl_mull(P17.x, inv_outputs[7]), ef_intrl_mull(P17.l, inv_outputs[7])};
	table[9] = (ec_point_laffine) {ef_intrl_mull(P19.x, inv_outputs[8]), ef_intrl_mull(P19.l, inv_outputs[8])};
	table[10] = (ec_point_laffine) {ef_intrl_mull(P21.x, inv_outputs[9]), ef_intrl_mull(P21.l, inv_outputs[9])};
	table[11] = (ec_point_laffine) {ef_intrl_mull(P23.x, inv_outputs[10]), ef_intrl_mull(P23.l, inv_outputs[10])};
	table[12] = (ec_point_laffine) {ef_intrl_mull(P25.x, inv_outputs[11]), ef_intrl_mull(P25.l, inv_outputs[11])};
	table[13] = (ec_point_laffine) {ef_intrl_mull(P27.x, inv_outputs[12]), ef_intrl_mull(P27.l, inv_outputs[12])};
	table[14] = (ec_point_laffine) {ef_intrl_mull(P29.x, inv_outputs[13]), ef_intrl_mull(P29.l, inv_outputs[13])};
	table[15] = (ec_point_laffine) {ef_intrl_mull(P31.x, inv_outputs[14]), ef_intrl_mull(P31.l, inv_outputs[14])};
}

#define MASK(B)                         (((uint64_t)1 << (B)) - 1)

void reg_rec(uint64x2_t k, uint64_t w, signed char rec[], uint64_t l) {
	uint64_t mask = MASK(w);

	for(int i = 0; i < l-1; i++) {
		int64_t k_i_temp = (k[0] & mask) - ((int64_t)1 << (w - 1)); // (k mod 16 only needs lower word)
		rec[i] = k_i_temp;
		k[0] -=k_i_temp;

		// First shift lower bits
		k[0] = k[0] >> (w-1);

		// Then higher bits with w-1 bits carry
		uint64_t carries = k[1] & mask;

		uint64_t shift_res;
		asm ("ROR %[res], %[input], %[s];"
		: [res] "=r" (shift_res)
		: [input] "r" (carries), [s] "r" (w-1)
		);

		k[0] = k[0] | shift_res;
		k[1] = k[1] >> (w-1);
	}
	rec[l-1] = k[0] & mask;
}

void ec_print_rec(signed char *rec, uint64_t l) {
  for(int i = l-1; i >= 0; i--) {
    printf(" %hhd ", rec[i]);
  }

  printf("\n");
}

/* GLS recoding in constant time, using new approach. */
ec_split_scalar ec_scalar_decomp(uint64x2x2_t k) {
	//c_i <= (q+1+|t|)/((q-1)^2 + t^2)*2^256 <= q+1+2sqrt(q)/(q+1)^2 * 2^256 <= 2^130-1
	//So the c_i will at most be 130 bits
    const uint64_t c1[4] = {0xE334618602D4CB44, 0x0000000000000004, 0x2, 0x0};
    const uint64_t c2[4] = {0x1CCB9E79FD2B34AC, 0xFFFFFFFFFFFFFFFB, 0x1, 0x0};
	//|a_i| <= q+1+|t| <= q+sqrt(q)+1 <= 2q-1
	//So the a_i will at most be 128 bits
    const uint64_t a1[2] = {0x38CD186180B532D2, 0x8000000000000001};
    const uint64_t a2[2] = {0xC732E79E7F4ACD2C, 0x7FFFFFFFFFFFFFFE};
    uint64_t v1[2][2], v2[2][2];
    uint64_t b1[2], b2[2], p1, p2, tmp1[4], tmp2[8];
	uint64_t zero = 0;
	ec_split_scalar res;

	//Step 1: Putting k into tmp1 for later multiplications
	tmp1[0] = k.val[0][0];
	tmp1[1] = k.val[0][1];
	tmp1[2] = k.val[1][0];
	tmp1[3] = k.val[1][1];

	//Step 2: Initializing basis
    /* v1 = (alpha1//2, alpha2//2), v2 = (alpha2//2, -alpha1//2) */
	//Need to remember that v2[1] is negative
    v2[1][0] = v1[0][0] = (a1[0] >> 1) | a1[1] << 63;
    v2[1][1] = v1[0][1] = (a1[1] >> 1);
    v2[0][0] = v1[1][0] = (a2[0] >> 1) | a2[1] << 63;
    v2[0][1] = v1[1][1] = (a2[1] >> 1);

    /* Step 3: Computing b1 = (c1*k) // 2^d, b2 = (c2*k) // 2^d for d = 256. */
    mult_u64_multiword(tmp2, c1, tmp1, 4);
    b1[0] = tmp2[4];
    b1[1] = tmp2[5];
    mult_u64_multiword(tmp2, c2, tmp1, 4);
    b2[0] = tmp2[4];
    b2[1] = tmp2[5];

	//Step 4: Compute parity corrections p1,p2 ahead of time
	//v2[0] is always even, so doesn't matter for parity.
	//b1*v1[0] is odd iff b1 is odd.
	p1 = (tmp1[0] + b1[0] + 1) & 0x1;
	ADDACC_128(p1, zero, b1[0], b1[1]);
	//v1[1] is always even, v2[1] odd.
	p2 = (b2[0] + 1) & 0x1;
	ADDACC_128(p2, zero, b2[0], b2[1]);
	// (u1, u2) = (v1, v2) as alpha1 % 4 != 0

	//Step 5: Compute k1
    mult_u64_multiword(tmp2, b1, v1[0], 2); //b1'*v1[0]
    mult_u64_multiword(tmp2+4, b2, v2[0], 2); //b2'*v2[0]
    ADDACC_256(tmp2[0], tmp2[1], tmp2[2], tmp2[3], tmp2[4], tmp2[5], tmp2[6], tmp2[7]); //b1'*v1[0] + b2'*v2[0]
	SUBACC_256(tmp2[4], tmp2[5], tmp2[6], tmp2[7], tmp1[0], tmp1[1], tmp1[2], tmp1[3]); //k - (b1'*v1[0] + b2'*v2[0])
	res.k1_sign = tmp1[3] != 0; //top half of tmp1 is either 0,0 or -1,-1 
	//Take two's complement if needed.
	res.k1[0] = tmp1[0] ^ (-res.k1_sign);
	res.k1[1] = tmp1[1] ^ (-res.k1_sign);
	ADDACC_128(res.k1_sign, zero, res.k1[0], res.k1[1]);

	//Step 6: Compute k2
	mult_u64_multiword(tmp1, b1, v1[1], 2); //b1'*v1[1]
	mult_u64_multiword(tmp2, b2, v2[1], 2); //-b2'*v2[1] 
	//b1'*v1[1] > 0 and b2'*v2[1] < 0.
	SUBACC_256(tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp2[0], tmp2[1], tmp2[2], tmp2[3]); //-b1'*v1[1] -b2'*v2[1]
	res.k2_sign = tmp2[3] != 0; //top half of tmp2 is either 0,0 or -1,-1 
	//Take two's complement if needed.
	res.k2[0] = tmp2[0] ^ (-res.k2_sign);
	res.k2[1] = tmp2[1] ^ (-res.k2_sign);
	ADDACC_128(res.k2_sign, zero, res.k2[0], res.k2[1]);

	return res;
}

void ec_precompute_w4_table2D_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ef_intrl_elem inv_inputs[16];
	ef_intrl_elem inv_outputs[16];
	
	ec_point_lproj tmp_table_proj[4];
	ec_point_lproj P2;
	ec_double_mixed_ptr(P, &P2);
	ec_laffine_to_lproj_ptr(P, &tmp_table_proj[0]);
	ec_add_mixed_unchecked_ptr(P, &P2, &tmp_table_proj[1]);
	ec_double_then_add_ptr(P, &P2, &tmp_table_proj[2]);
	ec_double_then_add_ptr(P, &tmp_table_proj[1], &tmp_table_proj[3]);
	inv_inputs[0] = tmp_table_proj[1].z;
	inv_inputs[1] = tmp_table_proj[2].z;
	inv_inputs[2] = tmp_table_proj[3].z;
	
	//Converting table to affine
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3, 0);
	ec_point_laffine tmp_table[4];
	tmp_table[0] = *P;
	for (int i = 1; i < 4; i++) {
		tmp_table[i] = (ec_point_laffine) {ef_intrl_mull(tmp_table_proj[i].x, inv_outputs[i-1]), ef_intrl_mull(tmp_table_proj[i].l, inv_outputs[i-1])};
	}
	
	//Computing table in projective coords
	ec_point_lproj table_proj[16];
	for (int j = 0; j < 4; j++) {
		ec_point_laffine jpoint = ec_endo_laffine(tmp_table[j]);
		for (int i = 0; i < 4; i++) {
			int index = 4*i + j;
			ec_add_laffine_unchecked_ptr(&tmp_table[i], &jpoint, &table_proj[index]);
			inv_inputs[index] = table_proj[index].z;
		}
	}
	
	//Converting table to affine
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 16, 0);
	for (int i = 0; i < 16; i++) {
		table[i] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}
}

void ec_precompute_w4_table2D_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ef_intrl_elem inv_inputs[16];
	ef_intrl_elem inv_outputs[16];
	ec_point_lproj tmp_table_proj[3];
	ec_point_lproj tmp;

	ec_triple_mixed_ptr(P, &tmp_table_proj[0]);
	ec_double_alt_ptr(&tmp_table_proj[0], &tmp);
	ec_add_sub_mixed_unchecked_ptr(P, &tmp, &tmp_table_proj[2], &tmp_table_proj[1]);
	ec_neg_mut(&tmp_table_proj[1]);

	inv_inputs[0] = tmp_table_proj[0].z;
	inv_inputs[1] = tmp_table_proj[1].z;
	inv_inputs[2] = tmp_table_proj[2].z;
	
	//Converting table to affine
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3, 0);
	ec_point_laffine tmp_table[4];
	tmp_table[0] = *P;
	for (int i = 1; i < 4; i++) {
		tmp_table[i] = (ec_point_laffine) {ef_intrl_mull(tmp_table_proj[i-1].x, inv_outputs[i-1]), ef_intrl_mull(tmp_table_proj[i-1].l, inv_outputs[i-1])};
	}
	
	//Computing table in projective coords
	ec_point_lproj table_proj[16];
	for(int i = 0; i < 4; i++) {
		int index = 5*i;
		ec_add_endo_laffine_unchecked_ptr(&tmp_table[i], &table_proj[index]);
		inv_inputs[index] = table_proj[index].z;
	}
	for(int j = 1; j < 4; j++) {
		ec_point_laffine jEndoP = ec_endo_laffine(tmp_table[j]);
		for(int i = 0; i < j; i++) {
			int indexij = 4*i+j;
			int indexji = 4*j+i;
			ec_add_sub_laffine_unchecked_ptr(&tmp_table[i], &jEndoP, &table_proj[indexij], &table_proj[indexji]);
			inv_inputs[indexij] = table_proj[indexij].z;
			inv_inputs[indexji] = table_proj[indexji].z;
		}
	}
	
	//Converting table to affine
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 16, 0);
	for (int i = 0; i < 16; i++) {
		table[i] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}

	for(int j = 1; j < 4; j++) {
		for(int i = 0; i < j; i++) {
			int indexji = 4*j+i;
			ec_cond_endo(&table[indexji], 1);
		}
	}
}

void ec_lookup_from_w4_table2D_ptr(ec_split_scalar *decomp, signed char rec_k1[], signed char rec_k2[], ec_point_laffine table[], int i, ec_point_laffine *next) {
	uint64x2x2_t lookup_data = ec_get_lookup_data_table2D(rec_k1[i], rec_k2[i], decomp, 4);
	lin_pass_w4_table2D(next, table, lookup_data.val[0][0]);
	
	ec_cond_endo(next, lookup_data.val[1][0]);
	//Cond neg
	next->l.val[0][0] ^= lookup_data.val[1][1];
}

void ec_scalarmull_single_endo_w4_table2D_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R) {
	uint64_t new_ptr, con = 1;
	int l = 43;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l);
	reg_rec(decomp.k2, 4, rec_k2, l);

	ec_point_laffine table[16];
	ec_precompute_w4_table2D_nonopt_ptr(P, table);
	
	ec_point_laffine next;
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, l-1, &next);
	
	ec_point_lproj Q, Qtmp1, Qtmp2;
	ec_laffine_to_lproj_ptr(&next, &Q);

	for(int i=l-2; i>=1; i--) {
		ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, i, &next);
		ec_double_alt_ptr(&Q, &Qtmp1);
		ec_double_alt_ptr(&Qtmp1, &Qtmp2);
		ec_double_then_add_ptr(&next, &Qtmp2, &Q);
	}
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &next);
	ec_double_alt_ptr(&Q, &Qtmp1);
	ec_double_alt_ptr(&Qtmp1, &Qtmp2);
	ec_double_alt_ptr(&Qtmp2, &Qtmp1);
	ec_add_mixed_ptr(&next, &Qtmp1, &Q);

	ec_lproj_to_laffine_ptr(&Q, R, 1);
}

void ec_lookup_from_w4_table2D_bulk_ptr(ec_split_scalar *decomp, signed char rec_k1[], signed char rec_k2[], ec_point_laffine table[], int i1, int i2, ec_point_laffine *P1, ec_point_laffine *P2) {
	//Code to find correct scalar values, as we can't redefine P as -P for negative k1 or k2. But we can, different points!!
	//Val comp is 2's complement conversion, but why the zero?
	uint64x2x2_t lookup_data1 = ec_get_lookup_data_table2D(rec_k1[i1], rec_k2[i1], decomp, 4);
	uint64x2x2_t lookup_data2 = ec_get_lookup_data_table2D(rec_k1[i2], rec_k2[i2], decomp, 4);

	lin_pass_w4_table2D_bulk(P1, P2, table, (uint64_t) lookup_data1.val[0][0], (uint64_t) lookup_data2.val[0][0]);

	ec_cond_endo(P1, lookup_data1.val[1][0]);
	ec_cond_endo(P2, lookup_data2.val[1][0]);

	//Cond neg
	P1->l.val[0][0] ^= lookup_data1.val[1][1];
	P2->l.val[0][0] ^= lookup_data2.val[1][1];
}

void ec_scalarmull_single_endo_w4_table2D_bulk_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R) {
	uint64_t new_ptr, con = 1;
	int l = 43;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l);
	reg_rec(decomp.k2, 4, rec_k2, l);

	ec_point_laffine table[16];
	ec_precompute_w4_table2D_ptr(P, table);
	
	ec_point_laffine P1, P2;
	ec_lookup_from_w4_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, l-1, l-2, &P1, &P2);

	ec_point_lproj Q, Qtmp1, Qtmp2;
	ec_double_mixed_ptr(&P1, &Qtmp1);
	ec_double_alt_ptr(&Qtmp1, &Qtmp2);
	ec_double_then_add_ptr(&P2, &Qtmp2, &Q);
	for(int i=l-3; i > 0; i -= 2) {
		ec_lookup_from_w4_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, i, i-1, &P1, &P2);
		ec_double_alt_ptr(&Q, &Qtmp1);
		ec_double_alt_ptr(&Qtmp1, &Qtmp2);
		ec_double_then_add_ptr(&P1, &Qtmp2, &Q);
		ec_double_alt_ptr(&Q, &Qtmp1);
		ec_double_alt_ptr(&Qtmp1, &Qtmp2);
		ec_double_then_add_ptr(&P2, &Qtmp2, &Q);
	}
	//Complete formula:
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &P1);
	ec_double_alt_ptr(&Q, &Qtmp1);
	ec_double_alt_ptr(&Qtmp1, &Qtmp2);
	ec_double_alt_ptr(&Qtmp2, &Qtmp1);
	ec_add_mixed_ptr(&P1, &Qtmp1, &Q);

	ec_lproj_to_laffine_ptr(&Q, R, 1);
}

void ec_precompute_w3_table2D_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ec_point_lproj P2_proj, P3_proj;
	ec_double_mixed_ptr(P, &P2_proj);
	ec_add_mixed_unchecked_ptr(P, &P2_proj, &P3_proj);
	ef_intrl_elem z3Inv = ef_intrl_inv(P3_proj.z, 0);
	ec_point_laffine P3 = (ec_point_laffine) {ef_intrl_mull(P3_proj.x, z3Inv), ef_intrl_mull(P3_proj.l, z3Inv)};

	ec_point_lproj table_proj[4];
	ec_add_endo_laffine_unchecked_ptr(P, &table_proj[0]);
	ec_add_endo_laffine_unchecked_ptr(&P3, &table_proj[3]);
	ec_cond_endo(&P3, 1);
	ec_add_sub_laffine_unchecked_ptr(P, &P3, &table_proj[1], &table_proj[2]);

	ef_intrl_elem inv_inputs[4];
	for(int i = 0; i < 4; i++) {
		inv_inputs[i] = table_proj[i].z;
	}
	ef_intrl_elem inv_outputs[4];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 4, 0);
	for (int i = 0; i < 4; i++) {
		table[i] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}

	ec_cond_endo(&table[2], 1);
}

void ec_precompute_w3_table2D_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ec_point_lproj table_proj[4];
	ec_point_lproj P2;//, P3, tmp;
	ec_add_endo_laffine_unchecked_ptr(P, &table_proj[0]);
	ec_double_mixed_ptr(P, &P2);
	
	ec_point_lproj psiP2 = ec_endo_lproj(P2);
	ec_add_unchecked_ptr(&table_proj[0], &psiP2, &table_proj[1]);
	ec_add_unchecked_ptr(&table_proj[0], &P2, &table_proj[2]);
	ec_add_unchecked_ptr(&table_proj[1], &P2, &table_proj[3]);

	//ec_add_mixed_unchecked_ptr(P, &P2, &P3);
	// ec_point_laffine psiP;
	// ec_endo_laffine_ptr(P, &psiP);
	// ec_add_sub_mixed_unchecked_ptr(&psiP, &P3, &table_proj[2], &tmp);
	// ec_neg_mut(&tmp);
	// table_proj[1] = ec_endo_lproj(tmp);
	// ec_add_unchecked_ptr(&P2, &table_proj[1], &table_proj[3]);

	ef_intrl_elem inv_inputs[4];
	for(int i = 0; i < 4; i++) {
		inv_inputs[i] = table_proj[i].z;
	}
	ef_intrl_elem inv_outputs[4];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 4, 0);
	for (int i = 0; i < 4; i++) {
		table[i] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}
}

void ec_scalarmull_single_endo_w3_table2D_bulk_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R) {
	uint64_t new_ptr, con = 1;
	int w = 3;
	int l = 64;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, w, rec_k1, l);
	reg_rec(decomp.k2, w, rec_k2, l);

	ec_point_laffine table[4];
	ec_precompute_w3_table2D_ptr(P, table);
	
	ec_point_laffine P1, P2;
	ec_lookup_from_w3_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, l-1, l-2, &P1, &P2);

	ec_point_lproj Q, Qtmp1, Qtmp2;
	ec_double_mixed_ptr(&P1, &Qtmp1);
	ec_double_then_add_ptr(&P2, &Qtmp1, &Q);

	for(int i=l-3; i > 0; i -= 2) {
		ec_lookup_from_w3_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, i, i-1, &P1, &P2);
		ec_double_alt_ptr(&Q, &Qtmp1);
		ec_double_then_add_ptr(&P1, &Qtmp1, &Q);
		ec_double_alt_ptr(&Q, &Qtmp1);
		if(i == 1) {
			ec_double_alt_ptr(&Qtmp1, &Qtmp2); 
			ec_add_mixed_ptr(&P2, &Qtmp2, &Q);
		} else {
			ec_double_then_add_ptr(&P2, &Qtmp1, &Q);
		}
	}

	ec_lproj_to_laffine_ptr(&Q, R, 1);
}