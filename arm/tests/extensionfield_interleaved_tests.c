#include "extensionfield_interleaved_tests.h"

#include "../common/extensionfield_interleaved.h"
#include "../common/utils.h"
#include <stdio.h>

void ef_intrl_interleave_test_example(test_ctr *ctr) {
	//Arrange
	ef_elem a = (ef_elem) {{{12, 34}, {56, 78}}};
	ef_intrl_elem expected = (ef_intrl_elem) {{{12, 56}, {34, 78}}}; 
	
	//Act
	ef_intrl_elem actual = ef_intrl_interleave(a);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "ef_intrl_interleave_test_example FAILED");
}

void ef_intrl_disentangle_test_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem a = (ef_intrl_elem) {{{789, 343}, {111, 222}}};
	ef_elem expected = (ef_elem) {{{789, 111}, {343, 222}}}; 
	
	//Act
	ef_elem actual = ef_intrl_disentangle(a);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "ef_intrl_disentangle_test_example FAILED");
}

void ef_intrl_disentangle_test_cancels_with_interleave(test_ctr *ctr) {
	//Arrange
	ef_elem a = (ef_elem) {{{12, 34}, {56, 78}}};
	
	//Act
	ef_elem res = ef_intrl_disentangle(ef_intrl_interleave(a));
	
	//Assert
	uint64_t correct = ef_intrl_equal(a, res);
	assert_true(correct, ctr, "ef_intrl_disentangle_test_cancels_with_interleave FAILED");
}

void ef_intrl_add_test_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem a = (ef_intrl_elem) {{{31, 2049}, {1, 15}}}; 
	ef_intrl_elem b = (ef_intrl_elem) {{{7, 1}, {1, 16}}}; 
	ef_intrl_elem expected = (ef_intrl_elem) {{{24, 2048}, {0, 31}}}; 
	
	//Act
	ef_intrl_elem actual = ef_intrl_add(a, b);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "ef_intrl_add_test_example FAILED");
}

void ef_intrl_red_test_example(test_ctr *ctr) {
	//Arrange
	//Shuffled representation of c1*u + c0 should be {c0[0], c1[0]}, {c0[1], c1[1]}, ...
	
	//(z^192 + z^191 + z^190)*u + (z^253 + z^128 + z^127 + 1)
	ef_intrl_elem_unred c = (ef_intrl_elem_unred) {{{1, 0}, {9223372036854775808U, 0}, {1, 13835058055282163712U}, {2305843009213693952, 1}}};
	//(z^127 + z^126 + z^65 + z^63 + z)*u + (z^127 + z^126 + z^125 + z^64 +z^62 + z +1) 
	ef_intrl_elem expected = (ef_intrl_elem) {{{4611686018427387907U, 9223372036854775810U}, {16140901064495857665U, 13835058055282163714U}}};
	
	//Act
	ef_intrl_elem actual = ef_intrl_red(c);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "ef_intrl_red_test_example FAILED");
}

void ef_intrl_square_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {1099511627777, 72057594037927936};
	poly64x2_t a1 = {1, 4611686022722355200};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^126 + z^96 + 1)u + (z^120 + z^40 + 1)
	poly64x2_t e0 = {2306405959167115266, 3459608938750738435};
	poly64x2_t e1 = {2305843009213693955, 3458764513820540931};
	ef_intrl_elem expected = ef_intrl_interleave(ef_create_elem(e0, e1)); 
	//(z^125 + z^124 + z^65 + z^64 + z^61 + z + 1)u + (z^125 + z^124 + z^113 + z^112 + z^80 + z^65 + z^64 + z^61 + z^49 + z)
	
	//Act
	ef_intrl_elem actual = ef_intrl_square(a);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "ef_intrl_square_test_example FAILED");
}

void ef_intrl_square_test_one_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_intrl_elem one = ef_intrl_interleave(ef_create_elem(one0, one1));
	
	//Act
	ef_intrl_elem result = ef_intrl_square(one);
	
	//Assert
	uint64_t correct = ef_intrl_equal(one, result);
	assert_true(correct, ctr, "ef_intrl_square_test_one_is_one FAILED");
}

void ef_intrl_square_test_zero_is_zero(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem zero = (ef_intrl_elem) {{{0, 0}, {0, 0}}};
	
	//Act
	ef_intrl_elem result = ef_intrl_square(zero);
	
	//Assert
	uint64_t correct = ef_intrl_equal(zero, result);
	assert_true(correct, ctr, "ef_intrl_square_test_zero_is_zero FAILED");
}

void ef_intrl_square_test_crosscheck_pmull(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {2251799813685249, 32768};
	poly64x2_t a1 = {17592186044416, 17179869186};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^98 + z^65+z^44)u + z^79 + z^51 + 1
	
	//Act
	ef_intrl_elem a_squared = ef_intrl_square(a);
	ef_intrl_elem a_squared_mull = ef_intrl_mull(a, a);
	
	//Assert
	uint64_t correct = ef_intrl_equal(a_squared, a_squared_mull);
	assert_true(correct, ctr, "ef_intrl_square_test_crosscheck_pmull FAILED");
}

void ef_intrl_square_test_freshmans_dream(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {2251799813685249, 32768};
	poly64x2_t a1 = {17592186044416, 17179869186};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^98 + z^65+z^44)u + z^79 + z^51 + 1
	poly64x2_t b0 = {1099511627777, 72057594037927936};
	poly64x2_t b1 = {1, 4611686022722355200};
	ef_intrl_elem b = ef_intrl_interleave(ef_create_elem(b0, b1)); //(z^126 + z^96 + 1)u + (z^120 + z^40 + 1)
	
	//Act
	ef_intrl_elem a_squared = ef_intrl_square(a);
	ef_intrl_elem b_squared = ef_intrl_square(b);
	ef_intrl_elem plus_after = ef_intrl_add(a_squared, b_squared);
	ef_intrl_elem a_plus_b = ef_intrl_add(a, b);
	ef_intrl_elem plus_before = ef_intrl_square(a_plus_b);
	
	//Assert
	uint64_t correct = ef_intrl_equal(plus_before, plus_after);
	assert_true(correct, ctr, "ef_intrl_square_test_freshmans_dream FAILED");
}

void ef_intrl_square_test_freshmans_dream_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		
		//Act
		ef_intrl_elem a_squared = ef_intrl_square(a);
		ef_intrl_elem b_squared = ef_intrl_square(b);
		ef_intrl_elem plus_after = ef_intrl_add(a_squared, b_squared);
		ef_intrl_elem a_plus_b = ef_intrl_add(a, b);
		ef_intrl_elem plus_before = ef_intrl_square(a_plus_b);
		
		//Assert
		correct &= ef_intrl_equal(plus_before, plus_after);
		if(!correct) {
			printf("a: ");
			ef_intrl_print_hex_nl(a);
			printf("b: ");
			ef_intrl_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "ef_intrl_square_test_freshmans_dream_rnd FAILED");
}

void ef_intrl_mull_A_test_crosscheck_ef_rnd(test_ctr *ctr) {
	//Arrange
	ef_elem a = ef_rand_elem();
	ef_intrl_elem a_intrl = ef_intrl_interleave(a);
	
	//Act
	ef_elem expected = ef_mull_A(a);
	ef_intrl_elem actual = ef_intrl_mull_A(a_intrl);
	
	//Assert
	uint64_t correct = ef_equal(expected, ef_intrl_disentangle(actual));
	assert_true(correct, ctr, "ef_intrl_mull_A_test_crosscheck_ef_rnd FAILED");
}

void ef_intrl_mull_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {9223372036854775809U, 9223372036854775808U}; 
	poly64x2_t a1 = {1, 1};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^64+1)u + (z^127 + z^63 + 1)
	poly64x2_t b0 = {4, 0};
	poly64x2_t b1 = {4, 4}; 
	ef_intrl_elem b = ef_intrl_interleave(ef_create_elem(b0, b1)); //(z^66 + z^2)u + (z^2)
	poly64x2_t e0 = {12, 4};
	poly64x2_t e1 = {8, 0};
	ef_intrl_elem expected = ef_intrl_interleave(ef_create_elem(e0, e1)); //(z^3)u + (z^66 + z^3 + z^2)
	
	//Act
	ef_intrl_elem actual = ef_intrl_mull(a, b);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "ef_intrl_mull_test_example FAILED");
}

void ef_intrl_mull_test_associative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {65792, 0};
	poly64x2_t a1 = {4294967296, 1};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^64 + z^32 )u + z^16 + z^8
	poly64x2_t b0 = {16, 0};
	poly64x2_t b1 = {9223372036854775808U, 0};
	ef_intrl_elem b = ef_intrl_interleave(ef_create_elem(b0, b1)); //(z^63)u + z^4
	poly64x2_t c0 = {36028797018963968, 70368744177664};
	poly64x2_t c1 = {0, 576460752303423488};
	ef_intrl_elem c = ef_intrl_interleave(ef_create_elem(c0, c1)); //(z^123)u + z^110 + z^55
	
	//Act
	ef_intrl_elem ab_first = ef_intrl_mull(ef_intrl_mull(a,b), c);
	ef_intrl_elem bc_first = ef_intrl_mull(a, ef_intrl_mull(b,c));
	
	//Assert
	uint64_t equal = ef_intrl_equal(ab_first, bc_first);
	assert_true(equal, ctr, "ef_intrl_mull_test_associative FAILED");
}

void ef_intrl_mull_test_associative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		ef_intrl_elem c = ef_intrl_rand_elem();
		
		//Act
		ef_intrl_elem ab_first = ef_intrl_mull(ef_intrl_mull(a,b), c);
		ef_intrl_elem bc_first = ef_intrl_mull(a, ef_intrl_mull(b,c));
		
		//Assert
		correct &= ef_intrl_equal(ab_first, bc_first);
		if(!correct) {
			printf("a: ");
			ef_intrl_print_hex_nl(a);
			printf("b: ");
			ef_intrl_print_hex_nl(b);
			printf("c: ");
			ef_intrl_print_hex_nl(c);
			break;
		}
	}
	assert_true(correct, ctr, "ef_intrl_mull_test_associative_rnd FAILED");
}

void ef_intrl_mull_test_commutative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {262144, 268435456};
	poly64x2_t a1 = {8589934592, 8192};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^77 + z^33)u + z^92 + z^18
	poly64x2_t b0 = {1025, 68719476736};
	poly64x2_t b1 = {2, 137438953473};
	ef_intrl_elem b = ef_intrl_interleave(ef_create_elem(b0, b1)); //(z^101 + z^64 + z)u + z^100 + z^10 + 1
	
	//Act
	ef_intrl_elem ab = ef_intrl_mull(a, b);
	ef_intrl_elem ba = ef_intrl_mull(b, a);
	
	//Assert
	uint64_t equal = ef_intrl_equal(ab, ba);
	assert_true(equal, ctr, "ef_intrl_mull_test_commutative FAILED");
}

void ef_intrl_mull_test_commutative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		
		//Act
		ef_intrl_elem ab = ef_intrl_mull(a, b);
		ef_intrl_elem ba = ef_intrl_mull(b, a);
		
		//Assert
		correct &= ef_intrl_equal(ab, ba);
		if(!correct) {
			printf("a: ");
			ef_intrl_print_hex_nl(a);
			printf("b: ");
			ef_intrl_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "ef_intrl_mull_test_commutative_rnd FAILED");
}

void ef_intrl_mull_test_one_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {40, 2305843009213693952};
	poly64x2_t a1 = {524288, 1048576};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^84 + z^19)u + (z^125 + z^5 + z^3)
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_intrl_elem one = ef_intrl_interleave(ef_create_elem(one0, one1));
	
	//Act
	ef_intrl_elem result = ef_intrl_mull(a, one);
	
	//Assert
	uint64_t correct = ef_intrl_equal(a, result);
	assert_true(correct, ctr, "ef_intrl_mull_test_one_is_identity FAILED");
}

void ef_intrl_mull_test_zero_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {4, 36028797018963968};
	poly64x2_t a1 = {17179869184, 34359738368};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^99 + z^34)u + (z^119 + z^2)
	poly64x2_t zero0 = {0, 0};
	ef_intrl_elem zero = ef_intrl_interleave(ef_create_elem(zero0, zero0));
	
	//Act
	ef_intrl_elem result = ef_intrl_mull(a, zero);
	
	//Assert
	uint64_t correct = ef_intrl_equal(zero, result);
	assert_true(correct, ctr, "ef_intrl_mull_test_zero_is_zero FAILED");
}

void ef_intrl_inv_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {4611686018427387904, 4611686018427387904};
	poly64x2_t a1 = {4611686018427387904, 4611686018427387904};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^126 + z^62)u + z^126 + z^62
	poly64x2_t e0 = {0, 0};
	poly64x2_t e1 = {2, 0};
	ef_intrl_elem expected = ef_intrl_interleave(ef_create_elem(e0, e1));
	
	//Act
	ef_intrl_elem actual = ef_intrl_inv(a);
	
	//Assert
	uint64_t correct = ef_intrl_equal(expected, actual);
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_example FAILED");
}

void ef_intrl_inv_test_inverse_of_one_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_intrl_elem one = ef_intrl_interleave(ef_create_elem(one0, one1));
	
	//Act
	ef_intrl_elem one_inv = ef_intrl_inv(one);
	
	//Assert
	uint64_t correct = ef_intrl_equal(one, one_inv);
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_inverse_of_one_is_one FAILED");
}

void ef_intrl_inv_zero_outputs_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t zero0 = {0, 0};
	ef_intrl_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_intrl_elem zero_inv = ef_intrl_inv(zero);
	
	//Assert
	uint64_t correct = ef_intrl_equal(zero, zero_inv);
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_zero_outputs_zero FAILED");
}

void ef_intrl_inv_test_inverse_of_inverse_is_original(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {8194, 144115188075855872};
	poly64x2_t a1 = {134217728, 16384};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^78 + z^27)u + (z^121 + z^13+z)
	
	//Act
	ef_intrl_elem a_inv = ef_intrl_inv(a);
	ef_intrl_elem a_inv_inv = ef_intrl_inv(a_inv);
	
	//Assert
	uint64_t correct = ef_intrl_equal(a, a_inv_inv);
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_inverse_of_inverse_is_original FAILED");
}

void ef_intrl_inv_test_inverse_of_inverse_is_original_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_intrl_elem a = ef_intrl_rand_elem();
		
		//Act
		ef_intrl_elem a_inv = ef_intrl_inv(a);
		ef_intrl_elem a_inv_inv = ef_intrl_inv(a_inv);
	
		//Assert
		correct &= ef_intrl_equal(a, a_inv_inv);
		if(!correct) {
			printf("a: ");
			ef_intrl_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_inverse_of_inverse_is_original_rnd FAILED");
}

void ef_intrl_inv_test_prod_of_inverses_is_inverse_of_prod(test_ctr *ctr) {
	//Assert
	poly64x2_t a0 = {40, 2305843009213693952};
	poly64x2_t a1 = {524288, 1048576};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^84 + z^19)u + (z^125 + z^5 + z^3)
	poly64x2_t b0 = {1025, 68719476736};
	poly64x2_t b1 = {2, 137438953473};
	ef_intrl_elem b = ef_intrl_interleave(ef_create_elem(b0, b1)); //(z^101 + z^64 + z)u + z^100 + z^10 + 1
	
	//Act
	ef_intrl_elem a_inv = ef_intrl_inv(a);
	ef_intrl_elem b_inv = ef_intrl_inv(b);
	ef_intrl_elem prod_of_inverses = ef_intrl_mull(a_inv, b_inv);
	ef_intrl_elem prod = ef_intrl_mull(a,b);
	ef_intrl_elem inv_of_prod = ef_intrl_inv(prod);
	
	//Assert
	uint64_t correct = ef_intrl_equal(prod_of_inverses, inv_of_prod);
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_prod_of_inverses_is_inverse_of_prod FAILED");
}
	
void ef_intrl_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		
		//Act
		ef_intrl_elem a_inv = ef_intrl_inv(a);
		ef_intrl_elem b_inv = ef_intrl_inv(b);
		ef_intrl_elem prod_of_inverses = ef_intrl_mull(a_inv, b_inv);
		ef_intrl_elem prod = ef_intrl_mull(a,b);
		ef_intrl_elem inv_of_prod = ef_intrl_inv(prod);
		
		//Assert
		correct = ef_intrl_equal(prod_of_inverses, inv_of_prod);
		if(!correct) {
			printf("a: ");
			ef_intrl_print_hex_nl(a);
			printf("b: ");
			ef_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_prod_of_inverses_is_inverse_of_prod_rnd FAILED");
}

void ef_intrl_inv_test_prod_with_inv_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {262144, 268435456};
	poly64x2_t a1 = {8589934592, 8192};
	ef_intrl_elem a = ef_intrl_interleave(ef_create_elem(a0, a1)); //(z^77 + z^33)u + z^92 + z^18
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_intrl_elem one = ef_intrl_interleave(ef_create_elem(one0, one1));
	
	//Act
	ef_intrl_elem a_inv = ef_intrl_inv(a);
	ef_intrl_elem a_times_inv = ef_intrl_mull(a, a_inv);
	
	ef_elem av = ef_intrl_disentangle(a);
	
	//Assert
	uint64_t correct = ef_intrl_equal(a_times_inv, one);
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_prod_with_inv_is_one FAILED");
}

void ef_intrl_inv_test_prod_with_inv_is_one_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_intrl_elem a = ef_intrl_rand_elem();
		poly64x2_t one0 = {1, 0};
		poly64x2_t one1 = {0, 0};
		ef_intrl_elem one = ef_intrl_interleave(ef_create_elem(one0, one1));
		
		//Act
		ef_intrl_elem a_inv = ef_intrl_inv(a);
		ef_intrl_elem a_times_inv = ef_intrl_mull(a, a_inv);
		
		//Assert
		correct = ef_intrl_equal(a_times_inv, one);
		if(!correct) {
			printf("a: ");
			ef_intrl_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_intrl_inv_test_prod_with_inv_is_one_rnd FAILED");
}

void ef_intrl_sim_inv_test_crosscheck_inv_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		uint64_t len = 3;
		ef_intrl_elem inputs[len];
		ef_intrl_elem outputs[len];
		inputs[0] = ef_intrl_rand_elem();
		inputs[1] = ef_intrl_rand_elem();
		inputs[2] = ef_intrl_rand_elem();
		
		//Act
		ef_intrl_sim_inv(inputs, outputs, len);
		
		//Assert
		correct = ef_intrl_equal(outputs[0], ef_intrl_inv(inputs[0])) 
			   && ef_intrl_equal(outputs[1], ef_intrl_inv(inputs[1])) 
			   && ef_intrl_equal(outputs[2], ef_intrl_inv(inputs[2]));
		if(!correct) {
			printf("Inputs[0]:\n");
			ef_intrl_print_hex_nl(inputs[0]);
			printf("Inputs[1]:\n");
			ef_intrl_print_hex_nl(inputs[1]);
			printf("Inputs[2]:\n");
			ef_intrl_print_hex_nl(inputs[2]);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_intrl_sim_inv_test_crosscheck_inv_rnd FAILED");
}

void extensionfield_interleaved_tests(test_ctr *ctr) {
	ef_intrl_interleave_test_example(ctr);
	ef_intrl_disentangle_test_example(ctr);
	ef_intrl_disentangle_test_cancels_with_interleave(ctr);
	
	ef_intrl_add_test_example(ctr);
	
	ef_intrl_red_test_example(ctr);
	
	ef_intrl_square_test_example(ctr);
	ef_intrl_square_test_one_is_one(ctr);
	ef_intrl_square_test_zero_is_zero(ctr);
	ef_intrl_square_test_crosscheck_pmull(ctr);
	ef_intrl_square_test_freshmans_dream(ctr);
	ef_intrl_square_test_freshmans_dream_rnd(ctr);
	
	ef_intrl_mull_A_test_crosscheck_ef_rnd(ctr);
	
	ef_intrl_mull_test_example(ctr);
	ef_intrl_mull_test_associative(ctr);
	ef_intrl_mull_test_associative_rnd(ctr);
	ef_intrl_mull_test_commutative(ctr);
	ef_intrl_mull_test_commutative_rnd(ctr);
	ef_intrl_mull_test_one_is_identity(ctr);
	ef_intrl_mull_test_zero_is_zero(ctr);
	
	ef_intrl_inv_test_example(ctr);
	ef_intrl_inv_test_inverse_of_one_is_one(ctr);
	ef_intrl_inv_zero_outputs_zero(ctr);
	ef_intrl_inv_test_inverse_of_inverse_is_original(ctr);
	ef_intrl_inv_test_inverse_of_inverse_is_original_rnd(ctr);
	ef_intrl_inv_test_prod_of_inverses_is_inverse_of_prod(ctr);
	ef_intrl_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(ctr);
	ef_intrl_inv_test_prod_with_inv_is_one(ctr);
	ef_intrl_inv_test_prod_with_inv_is_one_rnd(ctr);
	
	ef_intrl_sim_inv_test_crosscheck_inv_rnd(ctr);
}
