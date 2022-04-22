#include <stdio.h>

#include "basefield_tests.h"
#include "../common/basefield.h"
#include "../common/utils.h"

//========== Tests ===========================

void bf_create_elem_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t expected = { 86332548, 19205436754};
	
	//Act
	poly64x2_t actual = bf_create_elem(expected[0], expected[1]);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_create_elem_test_example FAILED");
}

void bf_add_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {4, 756};
	poly64x2_t b = {15, 2};
	poly64x2_t expected = {11, 758};
	
	//Act
	poly64x2_t actual = bf_add(a, b);
	
	//Assert
	uint64_t equal = equal_poly64x2(expected, actual); 
	assert_true(equal, ctr, "basefield: bf_add_test_example FAILED");
}

void bf_add_test_doubling_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {638540023, 896234759148761U};
	poly64x2_t zero = {0,0};
	
	//Act
	poly64x2_t result = bf_add(a, a);
	
	//Assert
	uint64_t is_zero = equal_poly64x2(result, zero);
	assert_true(is_zero, ctr, "basefield: bf_add_doubling_is_zero FAILED");
}

void bf_add_test_zero_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {129548752, 8236754276};
	poly64x2_t zero = {0,0};
	
	//Act
	poly64x2_t result = bf_add(a, zero);
	
	//Assert
	uint64_t is_identity = equal_poly64x2(result, a);
	assert_true(is_identity, ctr, "basefield: bf_add_zero_is_identity FAILED");
}

void bf_add_test_associative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {65192411524789, 444469420};
	poly64x2_t b = {3217681, 39129758209432};
	poly64x2_t c = {9437194765, 1111221111};
	
	//Act
	poly64x2_t abfirst = bf_add(bf_add(a, b), c);
	poly64x2_t bcfirst = bf_add(a, bf_add(b,c));
	
	//Assert
	uint64_t is_associative = equal_poly64x2(abfirst, bcfirst);
	assert_true(is_associative, ctr, "basefield: bf_add_test_associative FAILED");
}

void bf_add_test_commutative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {735761, 109136400836};
	poly64x2_t b = {5520459140631, 456059133};
	
	//Act
	poly64x2_t a_plus_b = bf_add(a, b);
	poly64x2_t b_plus_a = bf_add(b, a);
	
	//Assert
	uint64_t is_commutative = equal_poly64x2(a_plus_b, b_plus_a);
	assert_true(is_commutative, ctr, "basefield: bf_add_test_commutative FAILED");
}

void bf_pmull_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {1152921504606849024, 72057594037927936}; //z^120 + z^60 + z^11
	poly64x2_t b = {281474976714755, 70368744177664}; //z^110 + z^48 + z^12 + z + 1
	//z^230 + z^170 + z^168 + z^132 + z^120 + z^108 + z^72 + z^61 + z^60 + z^59 + z^23 + z^12 + z^11
	poly64x2_t expected0 = {4035225266132359168, 72075186223972608};
	poly64x2_t expected1 = {5497558138896, 274877906944};
	poly64x2x2_t expected = concat_bf_poly(expected0, expected1);
	
	//Act
	poly64x2x2_t actual = bf_pmull(a, b);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_pmull_test_example FAILED");
}

void bf_pmull_test_associative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {98314422310010283, 87159124330183471};
	poly64x2_t b = {3310980040412344311, 1120087461344321001};
	poly64x2_t c = {2276490184612864, 334652613944236613};
	
	//Act
	poly64x2_t ab = bf_red(bf_pmull(a,b));
	poly64x2_t bc = bf_red(bf_pmull(b,c));
	poly64x2_t ab_times_c = bf_red(bf_pmull(ab, c));
	poly64x2_t a_times_bc = bf_red(bf_pmull(a, bc));
	
	//Assert
	uint64_t equal = equal_poly64x2(ab_times_c, a_times_bc);
	assert_true(equal, ctr, "basefield: bf_pmull_test_associative FAILED");
}

void bf_pmull_test_associative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		poly64x2_t c = bf_rand_elem();
		
		//Act
		poly64x2_t ab = bf_red(bf_pmull(a,b));
		poly64x2_t bc = bf_red(bf_pmull(b,c));
		poly64x2_t ab_times_c = bf_red(bf_pmull(ab, c));
		poly64x2_t a_times_bc = bf_red(bf_pmull(a, bc));
		
		//Assert
		correct &= equal_poly64x2(ab_times_c, a_times_bc);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			printf("b: ");
			bf_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_pmull_test_associative_rnd FAILED");
}

void bf_pmull_test_commutative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {8, 549755817984}; //z^103 + z^76 + z^3
	poly64x2_t b = {2251799813685248, 16777216}; //z^89 + z^51 
	//Expected: z^192 + z^165 + z^154 + z^127 + z^92 + z^54;
	
	//Act
	poly64x2x2_t ab = bf_pmull(a,b);
	poly64x2x2_t ba = bf_pmull(b,a);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(ab, ba);
	assert_true(correct, ctr, "basefield: bf_pmull_test_commutative FAILED");
}

void bf_pmull_test_commutative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		
		//Act
		poly64x2x2_t ab = bf_pmull(a,b);
		poly64x2x2_t ba = bf_pmull(b,a);
		
		//Assert
		correct &= equal_poly64x2x2(ab, ba);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			printf("b: ");
			bf_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_pmull_test_commutative_rnd FAILED");
}

void bf_pmull_test_one_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {2251799813685248, 16777216}; //z^89 + z^51
	poly64x2_t one = {1,0};
	poly64x2_t zero = {0,0};
	poly64x2x2_t a_extended = concat_bf_poly(a, zero);
	
	//Act
	poly64x2x2_t res = bf_pmull(a, one);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(res, a_extended);
	assert_true(correct, ctr, "basefield: bf_pmull_test_one_is_identity FAILED");
}

void bf_pmull_test_zero_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {8, 549755817984}; //z^103 + z^76 + z^3
	poly64x2_t zero = {0,0};
	poly64x2x2_t zero_extended = concat_bf_poly(zero, zero);
	
	//Act
	poly64x2x2_t res = bf_pmull(a, zero);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(res, zero_extended);
	assert_true(correct, ctr, "basefield: bf_pmull_test_zero_is_zero FAILED");
}

void bf_psquare_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {8, 4611686018460942336}; //z^126 + z^89 + z^3
	
	poly64x2_t expected0 = {64, 0};
	poly64x2_t expected1 = {1125899906842624, 1152921504606846976};
	poly64x2x2_t expected = concat_bf_poly(expected0, expected1); //z^252 + z^178 + z^6
	
	//Act
	poly64x2x2_t actual = bf_psquare(a);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_psquare_test_example FAILED");
}

void bf_psquare_test_every_possible_term(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {18446744073709551615U, 9223372036854775807U}; //z^126 + z^125 + ... + z^2 + z + 1
	poly64x2_t e0 = {6148914691236517205, 6148914691236517205}; //z^126 + z^124 + ... + z^4 + z^2 + 1
	poly64x2_t e1 = {6148914691236517205, 1537228672809129301}; //z^252 + z^250 + ... + z^130 + z^128
	poly64x2x2_t expected = concat_bf_poly(e0, e1);
	
	//Act
	poly64x2x2_t actual = bf_psquare(a);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_psquare_test_every_possible_term FAILED");
}

void bf_psquare_test_one_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t one = {1, 0};
	poly64x2_t zero = {0, 0};
	poly64x2x2_t one_extended = concat_bf_poly(one, zero);
	
	//Act
	poly64x2x2_t result = bf_psquare(one);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(result, one_extended);
	assert_true(correct, ctr, "basefield: bf_psquare_test_one_is_one FAILED");
}

void bf_psquare_test_zero_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t zero = {0, 0};
	poly64x2x2_t zero_extended = concat_bf_poly(zero, zero);
	
	//Act
	poly64x2x2_t result = bf_psquare(zero);
	
	//Assert
	uint64_t correct = equal_poly64x2x2(result, zero_extended);
	assert_true(correct, ctr, "basefield: bf_psquare_test_zero_is_zero FAILED");
}

void bf_psquare_test_crosscheck_pmull(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {8796093087744, 36028797018963968}; //z^119 + z^43 + z^16
	
	//Act
	poly64x2x2_t result_psquare = bf_psquare(a);
	poly64x2x2_t result_pmull = bf_pmull(a, a);
	
	//Assert
	uint64_t equal = equal_poly64x2x2(result_psquare, result_pmull);
	assert_true(equal, ctr, "basefield: bf_psquare_test_crosscheck_pmull FAILED");
}

void bf_psquare_test_freshmans_dream(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {8796093087744, 36028797018963968}; //z^119 + z^43 + z^16
	poly64x2_t b = {4194308, 128}; //z^71 + z^22 + z^2
	
	//Act
	poly64x2x2_t apow2 = bf_psquare(a);
	poly64x2x2_t bpow2 = bf_psquare(b);
	poly64x2x2_t apow2_plus_bpow2 = concat_bf_poly(bf_add(apow2.val[0], bpow2.val[0]), bf_add(apow2.val[1], bpow2.val[1]));
	poly64x2x2_t aplusb_pow2 = bf_psquare(bf_add(a, b));
	
	//Assert
	uint64_t equal = equal_poly64x2x2(apow2_plus_bpow2, aplusb_pow2);
	assert_true(equal, ctr, "basefield: bf_psquare_test_freshmans_dream FAILED"); 
}

void bf_psquare_test_freshmans_dream_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		
		//Act
		poly64x2x2_t apow2 = bf_psquare(a);
		poly64x2x2_t bpow2 = bf_psquare(b);
		poly64x2x2_t apow2_plus_bpow2 = concat_bf_poly(bf_add(apow2.val[0], bpow2.val[0]), bf_add(apow2.val[1], bpow2.val[1]));
		poly64x2x2_t aplusb_pow2 = bf_psquare(bf_add(a, b));
		
		//Assert
		correct &= equal_poly64x2x2(apow2_plus_bpow2, aplusb_pow2);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			printf("b: ");
			bf_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_psquare_test_freshmans_dream_rnd FAILED");
}

void bf_red_test_example_twowords(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {36028797018963969, 13835058055282163712U};
	poly64x2_t zero = {0,0};
	poly64x2x2_t a = concat_bf_poly(a0, zero); //z^127 + z^126 + z^55 + 1
	poly64x2_t expected = {9259400833873739776U, 4611686018427387904}; //z^126 + z^63 + z^55
	
	//Act
	poly64x2_t actual = bf_red(a);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_red_test_example_twowords FAILED");
}

void bf_red_test_example_threewords(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {36028797018963969, 4611686018427387904U};
	poly64x2_t a1 = {2199023255554,0};
	poly64x2x2_t a = concat_bf_poly(a0, a1); //z^169 + z^129 + z^126 + z^55 + 1
	poly64x2_t expected = {36033195065475077, 4611688217450643458U}; 
	//z^126 + z^105 + z^65 + z^55 + z^42 + z^2 + 1
	
	//Act
	poly64x2_t actual = bf_red(a);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_red_test_example_threewords FAILED");
}

void bf_red_test_example_fourwords(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {1, 0};
	poly64x2_t a1 = {0, 1152921504606847040U};
	poly64x2x2_t a = concat_bf_poly(a0, a1); //z^252 + z^198 + 1
	poly64x2_t expected = {2305843009213694081U, 3458764513820541120U}; 
	//z^125 + z^124 + z^71 + z^70 + z^61 + z^7 + 1
	
	//Act
	poly64x2_t actual = bf_red(a);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_red_test_example_fourwords FAILED");
}

void bf_red_test_not_modifying_when_already_reduced(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {9223372036854775810U, 4611686018427387904U}; //z^126 + z^63 + z
	poly64x2_t zero = {0,0};
	poly64x2x2_t a_extended = concat_bf_poly(a, zero); 
	
	//Act
	poly64x2_t result = bf_red(a_extended);
	
	//Assert
	uint64_t unchanged = equal_poly64x2(result, a);
	assert_true(unchanged, ctr, "basefield: bf_red_test_not_modifying_when_already_reduced FAILED");
}

void bf_red_test_reduction_polynomial_reduces_to_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t f = {9223372036854775809U, 9223372036854775808U}; //z^127 + z^63 + 1
	poly64x2_t zero = {0,0};
	poly64x2x2_t f_extended = concat_bf_poly(f, zero);
	
	//Act
	poly64x2_t result = bf_red(f_extended);
	
	//Assert
	uint64_t iszero = equal_poly64x2(result, zero);
	assert_true(iszero, ctr, "basefield: bf_red_test_reduction_polynomial_reduces_to_zero FAILED");
}

void bf_red_psquare_test_example_threewords(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {4, 4611686018427387904};
	poly64x2_t a1 = {4611686018427387905, 0};
	poly64x2x2_t a = concat_bf_poly(a0, a1); // z^190 + z^128 + z^126 + z^2
	poly64x2_t expected = {9223372036854775814U, 1}; // z^64 + z^63 + z^2 + z
	
	//Act
	poly64x2_t actual = bf_red_psquare(a);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_red_psquare_test_example_threewords FAILED");
}

void bf_red_psquare_test_example_fourwords(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {16, 0};
	poly64x2_t a1 = {0, 1152921504606846976};
	poly64x2x2_t a = concat_bf_poly(a0, a1); //z^252 + z^4
	poly64x2_t expected = {2305843009213693968, 3458764513820540928}; //z^125 + z^124 + z^61 + z^4
	
	//Act
	poly64x2_t actual = bf_red_psquare(a);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_red_psquare_test_example_fourwords FAILED");
}

void bf_red_psquare_test_not_modifying_when_already_reduced(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {9223372036854775810U, 4611686018427387904U}; //z^126 + z^63 + z
	poly64x2_t zero = {0,0};
	poly64x2x2_t a_extended = concat_bf_poly(a, zero); 
	
	//Act
	poly64x2_t result = bf_red_psquare(a_extended);
	
	//Assert
	uint64_t unchanged = equal_poly64x2(result, a);
	assert_true(unchanged, ctr, "basefield: bf_red_psquare_test_not_modifying_when_already_reduced FAILED");
}

void bf_red_psquare_test_crosscheck_red(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {1, 1};
	poly64x2_t a1 = {4611686018427387905, 288230376151711745};
	poly64x2x2_t a = concat_bf_poly(a0, a1); // z^250 + z^192 + z^190 + z^128 + z^64 + 1
	
	//Act
	poly64x2_t a_red = bf_red(a);
	poly64x2_t a_red_psquare = bf_red_psquare(a);
	
	//Assert
	uint64_t equal = equal_poly64x2(a_red, a_red_psquare);
	assert_true(equal, ctr, "basefield: bf_red_psquare_test_crosscheck_red FAILED");
}

void bf_red_psquare_test_crosscheck_red_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		poly64x2x2_t a_squared = bf_psquare(a);
		
		//Act
		poly64x2_t a_squared_red = bf_red(a_squared);
		poly64x2_t a_squared_red_psquare = bf_red_psquare(a_squared);
		
		//Assert
		correct &= equal_poly64x2(a_squared_red, a_squared_red_psquare);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_red_psquare_test_crosscheck_red_rnd FAILED");
}

void bf_inv_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {0, 4611686018427387904U}; //z^126
	poly64x2_t expected = {3, 1}; //z^64 + z + 1
	
	//Act
	poly64x2_t actual = bf_inv(a);
	
	//Assert
	uint64_t correct = equal_poly64x2(expected, actual);
	assert_true(correct, ctr, "basefield: bf_inv_test_example FAILED");
}

void bf_inv_test_inverse_of_one_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t one = {1, 0};
	
	//Act
	poly64x2_t oneinv = bf_inv(one);
	
	//Assert
	uint64_t equal = equal_poly64x2(one, oneinv);
	assert_true(equal, ctr, "basefield: bf_inv_test_inverse_of_one_is_one FAILED");
}

void bf_inv_zero_outputs_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t zero = {0, 0};
	
	//Act
	poly64x2_t zeroinv = bf_inv(zero);
	
	//Assert
	uint64_t equal = equal_poly64x2(zero, zeroinv);
	assert_true(equal, ctr, "basefield: bf_inv_zero_outputs_zero FAILED");
}

void bf_inv_test_inverse_of_inverse_is_original(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {2233382993920, 68719477760}; //z^100 + z^74 + z^41 + z^35
	
	//Act
	poly64x2_t ainv = bf_inv(a);
	poly64x2_t ainvinv = bf_inv(ainv);
	
	//Assert
	uint64_t correct = equal_poly64x2(a, ainvinv);
	assert_true(correct, ctr, "basefield: bf_inv_test_inverse_of_inverse_is_original FAILED");
}

void bf_inv_test_inverse_of_inverse_is_original_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		
		//Act
		poly64x2_t ainv = bf_inv(a);
		poly64x2_t ainvinv = bf_inv(ainv);
	
		//Assert
		correct &= equal_poly64x2(a, ainvinv);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_inv_test_inverse_of_inverse_is_original_rnd FAILED");
}

void bf_inv_test_prod_of_inverses_is_inverse_of_prod(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {8, 549755817984}; //z^103 + z^76 + z^3
	poly64x2_t b = {2251799813685248, 16777216}; //z^89 + z^51 
	
	//Act
	poly64x2_t inv_of_prod = bf_inv(bf_red(bf_pmull(a,b)));
	poly64x2_t prod_of_inverses = bf_red(bf_pmull(bf_inv(a), bf_inv(b)));
	
	//Assert
	uint64_t correct = equal_poly64x2(inv_of_prod, prod_of_inverses);
	assert_true(correct, ctr, "basefield: bf_inv_test_prod_of_inverses_is_inverse_of_prod FAILED");
}

void bf_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		poly64x2_t zero = {0,0};
		
		if(equal_poly64x2(a, zero), equal_poly64x2(b, zero)) {
			continue;
		}
		
		//Act
		poly64x2_t inv_of_prod = bf_inv(bf_red(bf_pmull(a,b)));
		poly64x2_t prod_of_inverses = bf_red(bf_pmull(bf_inv(a), bf_inv(b)));
	
		//Assert
		correct &= equal_poly64x2(inv_of_prod, prod_of_inverses);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			printf("b: ");
			bf_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_inv_test_prod_of_inverses_is_inverse_of_prod_rnd FAILED");
}

void bf_inv_test_prod_with_inv_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {1125899940397057, 2305843077933172736U}; //z^125 + z^100 + z^75 + z^50 + z^25 + 1
	poly64x2_t one = {1, 0};
	
	//Act
	poly64x2_t ainv = bf_inv(a);
	poly64x2_t prod = bf_red(bf_pmull(a, ainv));
	
	//Assert
	uint64_t is_one = equal_poly64x2(prod, one);
	assert_true(is_one, ctr, "basefield: bf_inv_test_prod_with_inv_is_one FAILED");
}

void bf_inv_test_prod_with_inv_is_one_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		poly64x2_t a = bf_rand_elem();
		poly64x2_t one = {1, 0};
		poly64x2_t zero = {0,0};
		
		if(equal_poly64x2(a, zero)) {
			continue;
		}
	
		//Act
		poly64x2_t ainv = bf_inv(a);
		poly64x2_t prod = bf_red(bf_pmull(a, ainv));
	
		//Assert
		correct &= equal_poly64x2(prod, one);
		if(!correct) {
			printf("a: ");
			bf_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "basefield: bf_inv_test_prod_with_inv_is_one_rnd FAILED");
}

void basefield_tests(test_ctr *ctr) {
	bf_create_elem_test_example(ctr);
	
	bf_add_test_example(ctr);
	bf_add_test_doubling_is_zero(ctr);
	bf_add_test_zero_is_identity(ctr);
	bf_add_test_associative(ctr);
	bf_add_test_commutative(ctr);
	
	bf_pmull_test_example(ctr);
	bf_pmull_test_associative(ctr);
	bf_pmull_test_associative_rnd(ctr);
	bf_pmull_test_commutative(ctr);
	bf_pmull_test_commutative_rnd(ctr);
	bf_pmull_test_one_is_identity(ctr);
	bf_pmull_test_zero_is_zero(ctr);
	
	bf_psquare_test_example(ctr);
	bf_psquare_test_every_possible_term(ctr);
	bf_psquare_test_one_is_one(ctr);
	bf_psquare_test_zero_is_zero(ctr);
	bf_psquare_test_crosscheck_pmull(ctr);
	bf_psquare_test_freshmans_dream(ctr);
	bf_psquare_test_freshmans_dream_rnd(ctr);
	
	bf_red_test_example_twowords(ctr);
	bf_red_test_example_threewords(ctr);
	bf_red_test_example_fourwords(ctr);
	bf_red_test_not_modifying_when_already_reduced(ctr);
	bf_red_test_reduction_polynomial_reduces_to_zero(ctr);
	
	bf_red_psquare_test_example_threewords(ctr);
	bf_red_psquare_test_example_fourwords(ctr);
	bf_red_psquare_test_not_modifying_when_already_reduced(ctr);
	bf_red_psquare_test_crosscheck_red(ctr);
	bf_red_psquare_test_crosscheck_red_rnd(ctr);
	
	bf_inv_test_example(ctr);
	bf_inv_test_inverse_of_one_is_one(ctr);
	bf_inv_zero_outputs_zero(ctr);
	bf_inv_test_inverse_of_inverse_is_original(ctr);
	bf_inv_test_inverse_of_inverse_is_original_rnd(ctr);
	bf_inv_test_prod_of_inverses_is_inverse_of_prod(ctr);
	bf_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(ctr);
	bf_inv_test_prod_with_inv_is_one(ctr);
	bf_inv_test_prod_with_inv_is_one_rnd(ctr);
}
