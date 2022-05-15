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
	return ec_lproj_to_laffine(R);
}

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 65;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 3, rec_k1, l-1);
	reg_rec(decomp.k2, 3, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[2];
	precompute_w3(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
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
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));


	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 44;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l-1);
	reg_rec(decomp.k2, 4, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[4];
	precompute_w4(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
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
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w4(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 33;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l], rec_k2[l];

	reg_rec(decomp.k1, 5, rec_k1, l-1);
	reg_rec(decomp.k2, 5, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[8];
	precompute_w5(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1, P2;
	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
	//P1 = table[k1_val/2];
	//P2 = table[k2_val/2];

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
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
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

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));


	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

void ec_scalarmull_single_endo_w5_randaccess_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R) {
	uint64_t new_ptr, con = 1;
	int l = 33;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l], rec_k2[l];

	reg_rec(decomp.k1, 5, rec_k1, l-1);
	reg_rec(decomp.k2, 5, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[8];
	precompute_w5_ptr(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1, P2;
	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
	//P1 = table[k1_val/2];
	//P2 = table[k2_val/2];

	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q, Q1, Q2;
	ec_add_laffine_unchecked_ptr(&P1, &P2, &Q);

	for(int i=l-2; i>=0; i--) {
		ec_double_ptr(&Q, &Q1);
		ec_double_ptr(&Q1, &Q2);
		ec_double_ptr(&Q2, &Q1);

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		//P1 = table[k1_val/2];
		//P2 = table[k2_val/2];
		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		ec_double_then_addtwo_ptr(&P1, &P2, &Q1, &Q);
		//ec_print_hex_laffine(ec_lproj_to_laffine(Q));
	}

	P1 = ec_create_point_laffine(P->x, P->l);
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	P2 = ec_create_point_laffine(P->x, P->l);
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	ec_lproj_to_laffine_ptr(&Q, R);
}

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 27;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 6, rec_k1, l-1);
	reg_rec(decomp.k2, 6, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[16];
	precompute_w6(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
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
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
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

	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

void precompute_w3(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj Pl = ec_laffine_to_lproj(P);
	ec_point_lproj P3 = ec_double_then_add(P, Pl);

	ef_intrl_elem P3Z_inv = ef_intrl_inv(P3.z);

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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3);

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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 7);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
}

void precompute_w5_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ec_point_lproj P2;
	ec_double_mixed_ptr(P, &P2);
	ec_point_lproj P3;
	ec_add_mixed_unchecked_ptr(P, &P2, &P3);
	ec_point_lproj P4;
	ec_double_ptr(&P2, &P4);
	ec_point_lproj P5;
	ec_add_mixed_unchecked_ptr(P, &P4, &P5);
	ec_point_lproj P6; 
	ec_double_ptr(&P3, &P6);
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 7);

	table[0] = *P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 15);

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

  int i = 0;
  while (k[1] > 0 || i < l) {
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

    i++;
  }

	rec[i] = k[0] & mask;
}

void ec_print_rec(signed char *rec, uint64_t l) {
  for(int i = l-1; i >= 0; i--) {
    printf(" %hhd ", rec[i]);
  }

  printf("\n");
}

ec_split_scalar ec_scalar_decomp(uint64x2x2_t k) {
	ec_split_scalar result;
	uint64_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, sign;
	uint64_t zero = 0;

	// Step 1: {tmp1, tmp2, tmp3, tmp4} = b1 = k / 2^127 (where / is integer division)
	tmp0 = (k.val[1][0] << 1) | (k.val[0][1] >> 63);
	tmp1 = (k.val[1][1] << 1) | (k.val[1][0] >> 63);
	//tmp0 += ((k.val[0][1] >> 62) & 0x1); // round
	tmp2 = 0;
	tmp3 = 0;
	//printf("Step 1: %lu, %lu, %lu, %lu\n", tmp0, tmp1, tmp2, tmp3);

	//Step 2: b2 = k*trace / 2^254
	// Mult trace by each word in k.
	uint64x2_t tk0 = mult_u64((uint64_t) TRACE, k.val[0][0]);
	uint64x2_t tk1 = mult_u64((uint64_t) TRACE, k.val[0][1]);
	uint64x2_t tk2 = mult_u64((uint64_t) TRACE, k.val[1][0]);
	uint64x2_t tk3 = mult_u64((uint64_t) TRACE, k.val[1][1]);

	//Add previous [1] with next [0]
	//We are only interested in the last two words of the result, rest will be consumed by division by 2^254.
	asm volatile ("ADDS %[tmp4], %[tk01], %[tk10];"
		 "ADCS %[tmp4], %[tk11], %[tk20];"
		 "ADCS %[tmp4], %[tk21], %[tk30];"
		 "ADC %[tmp5], %[tk31], %[zero];"
		: [tmp4] "+r" (tmp4), [tmp5] "+r" (tmp5)
		: [tk01] "r" (tk0[1]), [tk10] "r" (tk1[0]), [tk11] "r" (tk1[1]), [tk20] "r" (tk2[0]), [tk21] "r" (tk2[1]), [tk30] "r" (tk3[0]), [tk31] "r" (tk3[1]), [zero] "r" (zero)
		);

	//Divide k*trace by 2^254
	uint64_t b2 = (tmp4 >> 62) | (tmp5 << 2);
	b2 += ((tmp4 >> 61) & 0x1); //Round
	//printf("Step 2: %lu\n", b2);

	//Step 3: b1*t
	uint64x2_t b10_times_t = mult_u64((uint64_t) TRACE, tmp0);
	uint64x2_t b11_times_t = mult_u64((uint64_t) TRACE, tmp1);
	uint64_t b1_times_t0, b1_times_t1, b1_times_t2;
	b1_times_t0 = b10_times_t[0];
	b1_times_t1 = b10_times_t[1];
	b1_times_t2 = b11_times_t[1];
	ADDACC_128(b11_times_t[0], zero, b1_times_t1, b1_times_t2);
	//printf("Step 3: %lu, %lu, %lu\n", b1_times_t0, b1_times_t1, b1_times_t2);

	//Step 4: b2*t
	uint64x2_t b2_times_t = mult_u64((uint64_t) TRACE, b2);
	//printf("Step 4: %lu, %lu\n", b2_times_t[0], b2_times_t[1]);

	//k1 computation

	//Step 5: {tmp4, tmp5, tmp6, tmp7} = b1*q (q = 2^127)
	tmp4 = 0;
	tmp5 = tmp0 << 63;
	tmp6 = (tmp0 >> 1) | (tmp1 << 63);
	tmp7 = tmp1 >> 1;
	//printf("Step 5: %lu, %lu, %lu, %lu\n", tmp4, tmp5, tmp6, tmp7);

	//Step 6: {tmp0, tmp1, tmp2, tmp3} = b1 + k
	ADDACC_256(k.val[0][0], k.val[0][1], k.val[1][0], k.val[1][1], tmp0, tmp1, tmp2, tmp3);
	//printf("%lu\n", k.val[0][0]);
	//printf("Step 6: %lu, %lu, %lu, %lu\n", tmp0, tmp1, tmp2, tmp3);

	//Step 7: {tmp4, tmp5, tmp6, tmp7} = b1*q + b2*t
	ADDACC_256(b2_times_t[0], b2_times_t[1], zero, zero, tmp4, tmp5, tmp6, tmp7);
	//printf("Step 7: %lu, %lu, %lu, %lu\n", tmp4, tmp5, tmp6, tmp7);

	//Step 8: Determine sign of k1 
	sign = (tmp3 < tmp7) | ((tmp3 == tmp7) & ((tmp2 < tmp6) | ((tmp2 == tmp6) & ((tmp1 < tmp5) | ((tmp1 == tmp5) & (tmp0 < tmp4))))));
	//printf("Step 8: %lu\n", sign);

	//Step 9: {tmp0, tmp1, tmp2, tmp3} = (b1 + k) - (b1*q + b2*t)
	SUBACC_256(tmp4, tmp5, tmp6, tmp7, tmp0, tmp1, tmp2, tmp3);
	//printf("Step 9: %lu, %lu, %lu, %lu\n", tmp0, tmp1, tmp2, tmp3);

	//Step 10: Now take two's complement if needed
	tmp0 ^= (zero - sign);
	tmp1 ^= (zero - sign);
	ADDACC_128(sign, zero, tmp0, tmp1);
	result.k1[0] = tmp0;
	result.k1[1] = tmp1;
	result.k1_sign = sign;
	//printf("Step 10: %lu, %lu, sign: %lu\n", tmp0, tmp1, sign);

	//k2 computation

	//Step 11: {tmp0, tmp1, tmp2} = b1*t + b2
	tmp0 = b1_times_t0;
	tmp1 = b1_times_t1;
	tmp2 = b1_times_t2;
	ADDACC_192(b2, zero, zero, tmp0, tmp1, tmp2);
	//printf("Step 11: %lu, %lu, %lu\n", tmp0, tmp1, tmp2);

	//Step 12: {tmp3, tmp4, tmp5} = b2*q (q = 2^127)
	tmp3 = 0;
	tmp4 = b2 << 63;
	tmp5 = b2 >> 1;
	//printf("Step 12: %lu, %lu, %lu\n", tmp3, tmp4, tmp5);

	//Step 13: k2 sign (0 for positive)
	sign = (tmp2 < tmp5) | ((tmp2 == tmp5) & ((tmp1 < tmp4) | ((tmp1 == tmp4) & (tmp0 < tmp3))));
	//printf("Step 13: %lu\n", sign);

	//Step 14: {tmp0, tmp1, tmp2} = b2*q - (b1*t + b2)
	SUBACC_192(tmp3, tmp4, tmp5, tmp0, tmp1, tmp2);
	//printf("Step 14: %lu, %lu, %lu\n", tmp0, tmp1, tmp2);

	//Step 15: Now take two's complement if needed.
	tmp0 ^= (zero - sign);
	tmp1 ^= (zero - sign);
	ADDACC_128(sign, zero, tmp0, tmp1);
	result.k2[0] = tmp0;
	result.k2[1] = tmp1;
	result.k2_sign = sign;
	//printf("Step 15: %lu, %lu sign: %lu\n", tmp0, tmp1, sign);

	return result;
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3);
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 16);
	for (int i = 0; i < 16; i++) {
		table[i] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}
}

void ec_precompute_w4_table2D_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3);
	ec_point_laffine tmp_table[4];
	tmp_table[0] = *P;
	for (int i = 1; i < 4; i++) {
		tmp_table[i] = (ec_point_laffine) {ef_intrl_mull(tmp_table_proj[i].x, inv_outputs[i-1]), ef_intrl_mull(tmp_table_proj[i].l, inv_outputs[i-1])};
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 16);
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
	int l = 44;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l-1);
	reg_rec(decomp.k2, 4, rec_k2, l-1);

	ec_point_laffine table[16];
	ec_precompute_w4_table2D_nonopt_ptr(P, table);
	
	ec_point_laffine next;
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, l-1, &next);
	//next = ec_lookup_from_w4_table2D(decomp, rec_k1, rec_k2, table, l-1);
	ec_point_lproj Q, Qtmp1, Qtmp2;
	ec_laffine_to_lproj_ptr(&next, &Q);
	// printf("Initial: \n");
	// ec_print_hex_laffine(ec_lproj_to_laffine(Q));

	for(int i=l-2; i>=0; i--) {
		ec_double_ptr(&Q, &Qtmp1);
		ec_double_ptr(&Qtmp1, &Qtmp2);
		ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, i, &next);
		//next = ec_lookup_from_w4_table2D(decomp, rec_k1, rec_k2, table, i);
		ec_double_then_add_ptr(&next, &Qtmp2, &Q);
		// ec_print_hex_laffine(ec_lproj_to_laffine(Q));
	}
	//Logic here with the xor is also strange
	ec_point_laffine P1 = ec_create_point_laffine(P->x, P->l);
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg;
	ec_add_mixed_unchecked_ptr(&P1, &Q, &Q_add_neg);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	ec_point_laffine P2 = ec_create_point_laffine(P->x, P->l);
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	ec_point_laffine psiP2 = ec_endo_laffine(P2);
	ec_add_mixed_unchecked_ptr(&psiP2, &Q, &Q_add_neg);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	ec_lproj_to_laffine_ptr(&Q, R); //Remember that this must use a constant time div
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
	int l = 44;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l-1);
	reg_rec(decomp.k2, 4, rec_k2, l-1);

	ec_point_laffine table[16];
	ec_precompute_w4_table2D_ptr(P, table);
	
	ec_point_laffine P1, P2;
	ec_lookup_from_w4_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, l-1, l-2, &P1, &P2);

	ec_point_lproj Q, Qtmp1, Qtmp2;
	ec_double_mixed_ptr(&P1, &Qtmp1);
	ec_double_ptr(&Qtmp1, &Qtmp2);
	ec_double_then_add_ptr(&P2, &Qtmp2, &Q);

	for(int i=l-3; i > 1; i -= 2) {
		ec_lookup_from_w4_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, i, i-1, &P1, &P2);
		ec_double_ptr(&Q, &Qtmp1);
		ec_double_ptr(&Qtmp1, &Qtmp2);
		ec_double_then_add_ptr(&P1, &Qtmp2, &Q);
		ec_double_ptr(&Q, &Qtmp1);
		ec_double_ptr(&Qtmp1, &Qtmp2);
		ec_double_then_add_ptr(&P2, &Qtmp2, &Q);
	}
	ec_lookup_from_w4_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, 1, 0, &P1, &P2);
	ec_double_ptr(&Q, &Qtmp1);
	ec_double_ptr(&Qtmp1, &Qtmp2);
	ec_double_then_add_ptr(&P1, &Qtmp2, &Q);
	ec_double_ptr(&Q, &Qtmp1);
	ec_double_ptr(&Qtmp1, &Qtmp2);
	ec_double_then_add_ptr(&P2, &Qtmp2, &Q); //Preparing for complete formula here

	//Logic here with the xor is also strange
	P1 = ec_create_point_laffine(P->x, P->l);
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg;
	ec_add_mixed_unchecked_ptr(&P1, &Q, &Q_add_neg);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	P2 = ec_create_point_laffine(P->x, P->l);
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	ec_point_laffine endoP2 = ec_endo_laffine(P2);
	ec_add_mixed_unchecked_ptr(&endoP2, &Q, &Q_add_neg);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	ec_lproj_to_laffine_ptr(&Q, R);
}

void ec_precompute_w3_table2D_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]) {
	ec_point_lproj P2_proj, P3_proj;
	ec_double_mixed_ptr(P, &P2_proj);
	ec_add_mixed_unchecked_ptr(P, &P2_proj, &P3_proj);
	ef_intrl_elem z3Inv = ef_intrl_inv(P3_proj.z);
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 4);
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
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 4);
	for (int i = 0; i < 4; i++) {
		table[i] = (ec_point_laffine) {ef_intrl_mull(table_proj[i].x, inv_outputs[i]), ef_intrl_mull(table_proj[i].l, inv_outputs[i])};
	}
}

void ec_scalarmull_single_endo_w3_table2D_bulk_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R) {
	uint64_t new_ptr, con = 1;
	int w = 3;
	int l = 65;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, w, rec_k1, l-1);
	reg_rec(decomp.k2, w, rec_k2, l-1);

	ec_point_laffine table[4];
	ec_precompute_w3_table2D_ptr(P, table);
	
	ec_point_laffine P1, P2;
	ec_lookup_from_w3_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, l-1, l-2, &P1, &P2);

	ec_point_lproj Q, Qtmp;
	ec_double_mixed_ptr(&P1, &Qtmp);
	ec_double_then_add_ptr(&P2, &Qtmp, &Q);
	for(int i=l-3; i > 0; i -= 2) {
		ec_lookup_from_w3_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, i, i-1, &P1, &P2);
		ec_double_ptr(&Q, &Qtmp);
		ec_double_then_add_ptr(&P1, &Qtmp, &Q);
		ec_double_ptr(&Q, &Qtmp);
		ec_double_then_add_ptr(&P2, &Qtmp, &Q);
	}
	ec_lookup_from_w3_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &P1);
	ec_double_ptr(&Q, &Qtmp);
	ec_double_then_add_ptr(&P1, &Qtmp, &Q); //Preparing for complete formula here

	//Logic here with the xor is also strange
	P1 = ec_create_point_laffine(P->x, P->l);
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg;
	ec_add_mixed_unchecked_ptr(&P1, &Q, &Q_add_neg);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	P2 = ec_create_point_laffine(P->x, P->l);
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	ec_point_laffine endoP2 = ec_endo_laffine(P2);
	ec_add_mixed_unchecked_ptr(&endoP2, &Q, &Q_add_neg);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	ec_lproj_to_laffine_ptr(&Q, R); //Remember this needs to be constant time
}