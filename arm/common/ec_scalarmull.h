#include <arm_neon.h>
#include "ec.h"
#include <stdio.h>

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

#include "utils.h"

typedef struct {
	uint64x2_t k1, k2;
	uint64_t k1_sign, k2_sign;
} ec_split_scalar;

extern void lin_pass_w3();

extern void lin_pass_w4();

extern void lin_pass_w5();

extern void lin_pass_w6();

extern void csel_asm();

extern void lin_pass_w3_table2D();

extern void lin_pass_w4_table2D();

extern void lin_pass_w3_table2D_bulk();

extern void lin_pass_w4_table2D_bulk();

ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l);

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k);

void ec_scalarmull_single_endo_w5_randaccess_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R);

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k);

void ec_scalarmull_single_endo_w4_table2D_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R);

void ec_scalarmull_single_endo_w3_table2D_bulk_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R);

void ec_scalarmull_single_endo_w4_table2D_bulk_ptr(ec_point_laffine *P, uint64x2x2_t k, ec_point_laffine *R);

void precompute_w3(ec_point_laffine P, ec_point_laffine* table);

void precompute_w4(ec_point_laffine P, ec_point_laffine* table);

void precompute_w5(ec_point_laffine P, ec_point_laffine* table);

void precompute_w5_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]);

void precompute_w5_ptr(ec_point_laffine *P, ec_point_laffine* table);

void precompute_w6(ec_point_laffine P, ec_point_laffine* table);

void ec_precompute_w4_table2D_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]);

void ec_precompute_w4_table2D_ptr(ec_point_laffine *P, ec_point_laffine table[]);

void ec_precompute_w3_table2D_nonopt_ptr(ec_point_laffine *P, ec_point_laffine table[]);

void ec_precompute_w3_table2D_ptr(ec_point_laffine *P, ec_point_laffine table[]);

void ec_lookup_from_w4_table2D_ptr(ec_split_scalar *decomp, signed char rec_k1[], signed char rec_k2[], ec_point_laffine table[], int i, ec_point_laffine *next);

void ec_lookup_from_w4_table2D_bulk_ptr(ec_split_scalar *decomp, signed char rec_k1[], signed char rec_k2[], ec_point_laffine table[], int i1, int i2, ec_point_laffine *P1, ec_point_laffine *P2);

void reg_rec(uint64x2_t k, uint64_t w, signed char* rec, uint64_t l);

void ec_print_rec(signed char *rec, uint64_t l);

ec_split_scalar ec_scalar_decomp(uint64x2x2_t k);

static inline uint64x2x2_t ec_get_lookup_data_table2D(signed char k1_digit, signed char k2_digit, ec_split_scalar *decomp, int ncols) {
	//Code to find correct scalar values, as we can't redefine P as -P for negative k1 or k2. But we can, different points!!
	//Val comp is 2's complement conversion, but why the zero?
	uint64_t zero = 0;
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign; //why zero
	uint64_t k1_sign = k1_digit_sign^decomp->k1_sign;

	uint64_t k2_digit_sign = ((unsigned char)k2_digit >> 7);
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp->k2_sign;
	
	uint64_t sign_xor = k1_sign ^ k2_sign;
	uint64_t normal_index = ncols*(k1_val/2)+(k2_val/2); //redundancy to round down
	uint64_t swap_index = ncols*(k2_val/2)+(k1_val/2);
	uint64_t index = sign_xor * swap_index + (1-sign_xor) * normal_index; 
	
	uint64x2x2_t res;
	res.val[0][0] = index;
	res.val[1][0] = sign_xor; //endo cond
	res.val[1][1] = k2_sign; //negation cond
	return res;
}

static inline void ec_cond_endo(ec_point_laffine *P, uint64_t cond) {
	P->x.val[0][0] ^= P->x.val[0][1]*cond;
	P->x.val[1][0] ^= P->x.val[1][1]*cond;
	P->l.val[0][0] ^= P->l.val[0][1]*cond;
	P->l.val[1][0] ^= P->l.val[1][1]*cond;
	P->l.val[0][1] ^= cond; 
}

static inline void ec_lookup_from_w3_table2D_ptr(ec_split_scalar *decomp, signed char rec_k1[], signed char rec_k2[], ec_point_laffine table[], int i, ec_point_laffine *P1) {
	uint64x2x2_t lookup_data = ec_get_lookup_data_table2D(rec_k1[i], rec_k2[i], decomp, 2);

	lin_pass_w3_table2D(P1, table, lookup_data.val[0][0]);

	ec_cond_endo(P1, lookup_data.val[1][0]);

	//Cond neg
	P1->l.val[0][0] ^= lookup_data.val[1][1];
}

static inline void ec_lookup_from_w3_table2D_bulk_ptr(ec_split_scalar *decomp, signed char rec_k1[], signed char rec_k2[], ec_point_laffine table[], int i1, int i2, ec_point_laffine *P1, ec_point_laffine *P2) {
	uint64x2x2_t lookup_data1 = ec_get_lookup_data_table2D(rec_k1[i1], rec_k2[i1], decomp, 2);
	uint64x2x2_t lookup_data2 = ec_get_lookup_data_table2D(rec_k1[i2], rec_k2[i2], decomp, 2);

	lin_pass_w3_table2D_bulk(P1, P2, table, lookup_data1.val[0][0], lookup_data2.val[0][0]);

	ec_cond_endo(P1, lookup_data1.val[1][0]);
	ec_cond_endo(P2, lookup_data2.val[1][0]);

	//Cond neg
	P1->l.val[0][0] ^= lookup_data1.val[1][1];
	P2->l.val[0][0] ^= lookup_data2.val[1][1];
}

#define all1s 18446744073709551615U //2^64 - 1

extern void lin_pass_w4_table2D_bulk_neon();

#endif
