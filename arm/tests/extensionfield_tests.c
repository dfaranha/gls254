#include <stdio.h>

#include "extensionfield_tests.h"
#include "../common/extensionfield.h"
#include "../common/utils.h"

void ef_create_elem_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {16785442, 721476};
	poly64x2_t a1 = {78099554548664, 6547959942615};
	
	//Act
	ef_elem a = ef_create_elem(a0, a1);
	
	//Assert
	uint64_t correct = equal_poly64x2(a.val[0], a0) & equal_poly64x2(a.val[1], a1);
	assert_true(correct, ctr, "extensionfield: ef_create_elem_test_example FAILED");
}

void ef_equal_test_equal(test_ctr *ctr) {
	//Arrange
	ef_elem a;
	poly64x2_t a0 = {6749534, 697497494916};
	poly64x2_t a1 = {7681514, 14870085};
	a.val[0] = a0;
	a.val[1] = a1;
	ef_elem b;
	poly64x2_t b0 = {6749534, 697497494916};
	poly64x2_t b1 = {7681514, 14870085};
	b.val[0] = b0;
	b.val[1] = b1;
	
	//Act
	uint64_t equal = ef_equal(a, b);
	
	//Assert
	assert_true(equal, ctr, "utils: ef_equal_test_equal FAILED");
}

void ef_equal_test_notequal(test_ctr *ctr) {
	//Arrange
	ef_elem a;
	poly64x2_t a0 = {6749534, 697497494916};
	poly64x2_t a1 = {7681514, 14870085};
	a.val[0] = a0;
	a.val[1] = a1;
	ef_elem b;
	poly64x2_t b0 = {6749534, 697497394916};
	poly64x2_t b1 = {7681514, 14870085};
	b.val[0] = b0;
	b.val[1] = b1;
	
	//Act
	uint64_t equal = ef_equal(a, b);
	
	//Assert
	assert_false(equal, ctr, "utils: ef_equal_test_notequal FAILED");
}

void ef_add_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {31, 7};
	poly64x2_t a1 = {0, 64};
	ef_elem a = ef_create_elem(a0, a1);
	
	poly64x2_t b0 = {9, 13};
	poly64x2_t b1 = {1, 64};
	ef_elem b = ef_create_elem(b0, b1);
	
	poly64x2_t e0 = {22, 10};
	poly64x2_t e1 = {1, 0};
	ef_elem expected = ef_create_elem(e0, e1);
	
	//Act
	ef_elem actual = ef_add(a, b);
	
	//Assert
	uint64_t correct = ef_equal(expected, actual);
	assert_true(correct, ctr, "extensionfield: ef_add_test_example FAILED");
}

void ef_add_test_doubling_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {6751410153, 112308523};
	poly64x2_t a1 = {998742121496553186, 718656530};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t zero0 = {0,0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem actual = ef_add(a, a);
	
	//Assert
	uint64_t is_zero = ef_equal(zero, actual);
	assert_true(is_zero, ctr, "extensionfield: ef_add_test_doubling_is_zero FAILED");
}

void ef_add_test_zero_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {233207851, 6715398730024};
	poly64x2_t a1 = {67651597459761, 12435676};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t zero0 = {0,0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem result = ef_add(a, zero);
	
	//Assert
	uint64_t correct = ef_equal(a, result);
	assert_true(correct, ctr, "extensionfield: ef_add_test_zero_is_identity FAILED");
}

void ef_add_test_associative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {7651232241, 2233445678};
	poly64x2_t a1 = {101, 59874};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t b0 = {149073,8569319135};
	poly64x2_t b1 = {2287569446512,433318594};
	ef_elem b = ef_create_elem(b0, b1);
	poly64x2_t c0 = {1786054612, 786156412};
	poly64x2_t c1 = {54364258769123, 3521913758};
	ef_elem c = ef_create_elem(c0, c1);
	
	//Act
	ef_elem a_plus_b_first = ef_add(ef_add(a, b), c);
	ef_elem b_plus_c_first = ef_add(a, ef_add(b, c));
	
	//Assert
	uint64_t equal = ef_equal(a_plus_b_first, b_plus_c_first);
	assert_true(equal, ctr, "extensionfield: ef_add_test_associative FAILED");
}

void ef_add_test_commutative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {9785313565, 7656548341100};
	poly64x2_t a1 = {94271531642, 2176487654};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t b0 = {78695976431414, 5677568658888};
	poly64x2_t b1 = {77012835631, 2176456423658};
	ef_elem b = ef_create_elem(b0, b1);
	
	//Act
	ef_elem a_plus_b = ef_add(a, b);
	ef_elem b_plus_a = ef_add(b, a);
	
	//Assert
	uint64_t equal = ef_equal(a_plus_b, b_plus_a);
	assert_true(equal, ctr, "extensionfield: ef_add_test_commutative FAILED");
}
	
void ef_mull_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {9223372036854775809U, 9223372036854775808U}; 
	poly64x2_t a1 = {1, 1};
	ef_elem a = ef_create_elem(a0, a1); //(z^64+1)u + (z^127 + z^63 + 1)
	poly64x2_t b0 = {4, 0};
	poly64x2_t b1 = {4, 4}; 
	ef_elem b = ef_create_elem(b0, b1); //(z^66 + z^2)u + (z^2)
	poly64x2_t e0 = {12, 4};
	poly64x2_t e1 = {8, 0};
	ef_elem expected = ef_create_elem(e0, e1); //(z^3)u + (z^66 + z^3 + z^2)
	
	//Act
	ef_elem actual = ef_mull(a, b);
	
	//Assert
	uint64_t correct = ef_equal(expected, actual);
	assert_true(correct, ctr, "extensionfield: ef_mull_test_example FAILED");
}

void ef_mull_test_associative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {65792, 0};
	poly64x2_t a1 = {4294967296, 1};
	ef_elem a = ef_create_elem(a0, a1); //(z^64 + z^32 )u + z^16 + z^8
	poly64x2_t b0 = {16, 0};
	poly64x2_t b1 = {9223372036854775808U, 0};
	ef_elem b = ef_create_elem(b0, b1); //(z^63)u + z^4
	poly64x2_t c0 = {36028797018963968, 70368744177664};
	poly64x2_t c1 = {0, 576460752303423488};
	ef_elem c = ef_create_elem(c0, c1); //(z^123)u + z^110 + z^55
	
	//Act
	ef_elem ab_first = ef_mull(ef_mull(a,b), c);
	ef_elem bc_first = ef_mull(a, ef_mull(b,c));
	
	//Assert
	uint64_t equal = ef_equal(ab_first, bc_first);
	assert_true(equal, ctr, "extensionfield: ef_mull_test_associative FAILED");
}

void ef_mull_test_associative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_elem a = ef_rand_elem();
		ef_elem b = ef_rand_elem();
		ef_elem c = ef_rand_elem();
		
		//Act
		ef_elem ab_first = ef_mull(ef_mull(a,b), c);
		ef_elem bc_first = ef_mull(a, ef_mull(b,c));
		
		//Assert
		correct &= ef_equal(ab_first, bc_first);
		if(!correct) {
			printf("a: ");
			ef_print_hex_nl(a);
			printf("b: ");
			ef_print_hex_nl(b);
			printf("c: ");
			ef_print_hex_nl(c);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_mull_test_associative_rnd FAILED");
}

void ef_mull_test_commutative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {262144, 268435456};
	poly64x2_t a1 = {8589934592, 8192};
	ef_elem a = ef_create_elem(a0, a1); //(z^77 + z^33)u + z^92 + z^18
	poly64x2_t b0 = {1025, 68719476736};
	poly64x2_t b1 = {2, 137438953473};
	ef_elem b = ef_create_elem(b0, b1); //(z^101 + z^64 + z)u + z^100 + z^10 + 1
	
	//Act
	ef_elem ab = ef_mull(a, b);
	ef_elem ba = ef_mull(b, a);
	
	//Assert
	uint64_t equal = ef_equal(ab, ba);
	assert_true(equal, ctr, "extensionfield: ef_mull_test_commutative FAILED");
}

void ef_mull_test_commutative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_elem a = ef_rand_elem();
		ef_elem b = ef_rand_elem();
		
		//Act
		ef_elem ab = ef_mull(a, b);
		ef_elem ba = ef_mull(b, a);
		
		//Assert
		correct &= ef_equal(ab, ba);
		if(!correct) {
			printf("a: ");
			ef_print_hex_nl(a);
			printf("b: ");
			ef_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_mull_test_commutative_rnd FAILED");
}

void ef_mull_test_one_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {40, 2305843009213693952};
	poly64x2_t a1 = {524288, 1048576};
	ef_elem a = ef_create_elem(a0, a1); //(z^84 + z^19)u + (z^125 + z^5 + z^3)
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_elem one = ef_create_elem(one0, one1);
	
	//Act
	ef_elem result = ef_mull(a, one);
	
	//Assert
	uint64_t correct = ef_equal(a, result);
	assert_true(correct, ctr, "extensionfield: ef_mull_test_one_is_identity FAILED");
}

void ef_mull_test_zero_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {4, 36028797018963968};
	poly64x2_t a1 = {17179869184, 34359738368};
	ef_elem a = ef_create_elem(a0, a1); //(z^99 + z^34)u + (z^119 + z^2)
	poly64x2_t zero0 = {0, 0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem result = ef_mull(a, zero);
	
	//Assert
	uint64_t correct = ef_equal(zero, result);
	assert_true(correct, ctr, "extensionfield: ef_mull_test_zero_is_zero FAILED");
}

void ef_mull_B_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {137438953472, 68719476736}; //x^100 + x^37
	poly64x2_t a1 = {274877906944, 13835058192721117184U}; //x^127 + x^126 + x^101 + x^38
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t e0 = {137438953472, 9223372105574252545U}; //x^127 + x^100 + x^64 + x^37
	poly64x2_t e1 = {275079233538, 13835058192821780483U}; // x^127 + x^126 + x^101 + x^90 + x^89 + x^65 + x^64 + x^38 + x^27 + x^26 + x
	ef_elem expected = ef_create_elem(e0, e1);
	
	//Act
	ef_elem actual = ef_mull_B(a);
	
	//Assert
	uint64_t correct = ef_equal(expected, actual);
	assert_true(correct, ctr, "ef_mull_B_test_example FAILED");
}
	
void ef_square_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {1099511627777, 72057594037927936};
	poly64x2_t a1 = {1, 4611686022722355200};
	ef_elem a = ef_create_elem(a0, a1); //(z^126 + z^96 + 1)u + (z^120 + z^40 + 1)
	poly64x2_t e0 = {2306405959167115266, 3459608938750738435};
	poly64x2_t e1 = {2305843009213693955, 3458764513820540931};
	ef_elem expected = ef_create_elem(e0, e1); //(z^125 + z^124 + z^65 + z^64 + z^61 + z + 1)u + (z^125 + z^124 + z^113 + z^112 + z^80 + z^65 + z^64 + z^61 + z^49 + z)
	
	//Act
	ef_elem actual = ef_square(a);
	
	//Assert
	uint64_t correct = ef_equal(expected, actual);
	assert_true(correct, ctr, "extensionfield: ef_square_test_example FAILED");
}

void ef_square_test_one_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_elem one = ef_create_elem(one0, one1);
	
	//Act
	ef_elem result = ef_square(one);
	
	//Assert
	uint64_t correct = ef_equal(one, result);
	assert_true(correct, ctr, "extensionfield: ef_square_test_one_is_one FAILED");
}

void ef_square_test_zero_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t zero0 = {0, 0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem result = ef_square(zero);
	
	//Assert
	uint64_t correct = ef_equal(zero, result);
	assert_true(correct, ctr, "extensionfield: ef_square_test_zero_is_zero FAILED");
}

void ef_square_test_crosscheck_pmull(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {2251799813685249, 32768};
	poly64x2_t a1 = {17592186044416, 17179869186};
	ef_elem a = ef_create_elem(a0, a1); //(z^98 + z^65+z^44)u + z^79 + z^51 + 1
	
	//Act
	ef_elem a_squared = ef_square(a);
	ef_elem a_squared_mull = ef_mull(a, a);
	
	//Assert
	uint64_t correct = ef_equal(a_squared, a_squared_mull);
	assert_true(correct, ctr, "extensionfield: ef_square_test_crosscheck_pmull FAILED");
}

void ef_square_test_freshmans_dream(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {2251799813685249, 32768};
	poly64x2_t a1 = {17592186044416, 17179869186};
	ef_elem a = ef_create_elem(a0, a1); //(z^98 + z^65+z^44)u + z^79 + z^51 + 1
	poly64x2_t b0 = {1099511627777, 72057594037927936};
	poly64x2_t b1 = {1, 4611686022722355200};
	ef_elem b = ef_create_elem(b0, b1); //(z^126 + z^96 + 1)u + (z^120 + z^40 + 1)
	
	//Act
	ef_elem a_squared = ef_square(a);
	ef_elem b_squared = ef_square(b);
	ef_elem plus_after = ef_add(a_squared, b_squared);
	ef_elem a_plus_b = ef_add(a, b);
	ef_elem plus_before = ef_square(a_plus_b);
	
	//Assert
	uint64_t correct = ef_equal(plus_before, plus_after);
	assert_true(correct, ctr, "extensionfield: ef_square_test_freshmans_dream FAILED");
}

void ef_square_test_freshmans_dream_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_elem a = ef_rand_elem();
		ef_elem b = ef_rand_elem();
		
		//Act
		ef_elem a_squared = ef_square(a);
		ef_elem b_squared = ef_square(b);
		ef_elem plus_after = ef_add(a_squared, b_squared);
		ef_elem a_plus_b = ef_add(a, b);
		ef_elem plus_before = ef_square(a_plus_b);
		
		//Assert
		correct &= ef_equal(plus_before, plus_after);
		if(!correct) {
			printf("a: ");
			ef_print_hex_nl(a);
			printf("b: ");
			ef_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_square_test_freshmans_dream_rnd FAILED");
}

void ef_inv_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {4611686018427387904, 4611686018427387904};
	poly64x2_t a1 = {4611686018427387904, 4611686018427387904};
	ef_elem a = ef_create_elem(a0, a1); //(z^126 + z^62)u + z^126 + z^62
	poly64x2_t e0 = {0, 0};
	poly64x2_t e1 = {2, 0};
	ef_elem expected = ef_create_elem(e0, e1);
	
	//Act
	ef_elem actual = ef_inv(a);
	
	//Assert
	uint64_t correct = ef_equal(expected, actual);
	assert_true(correct, ctr, "extensionfield: ef_inv_test_example FAILED");
}

void ef_inv_test_inverse_of_one_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_elem one = ef_create_elem(one0, one1);
	
	//Act
	ef_elem one_inv = ef_inv(one);
	
	//Assert
	uint64_t correct = ef_equal(one, one_inv);
	assert_true(correct, ctr, "extensionfield: ef_inv_test_inverse_of_one_is_one FAILED");
}

void ef_inv_zero_outputs_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t zero0 = {0, 0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem zero_inv = ef_inv(zero);
	
	//Assert
	uint64_t correct = ef_equal(zero, zero_inv);
	assert_true(correct, ctr, "extensionfield: ef_inv_zero_outputs_zero FAILED");
}

void ef_inv_test_inverse_of_inverse_is_original(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {8194, 144115188075855872};
	poly64x2_t a1 = {134217728, 16384};
	ef_elem a = ef_create_elem(a0, a1); //(z^78 + z^27)u + (z^121 + z^13+z)
	
	//Act
	ef_elem a_inv = ef_inv(a);
	ef_elem a_inv_inv = ef_inv(a_inv);
	
	//Assert
	uint64_t correct = ef_equal(a, a_inv_inv);
	assert_true(correct, ctr, "extensionfield: ef_inv_test_inverse_of_inverse_is_original FAILED");
}

void ef_inv_test_inverse_of_inverse_is_original_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_elem a = ef_rand_elem();
		
		//Act
		ef_elem a_inv = ef_inv(a);
		ef_elem a_inv_inv = ef_inv(a_inv);
	
		//Assert
		correct &= ef_equal(a, a_inv_inv);
		if(!correct) {
			printf("a: ");
			ef_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_inv_test_inverse_of_inverse_is_original_rnd FAILED");
}

void ef_inv_test_prod_of_inverses_is_inverse_of_prod(test_ctr *ctr) {
	//Assert
	poly64x2_t a0 = {40, 2305843009213693952};
	poly64x2_t a1 = {524288, 1048576};
	ef_elem a = ef_create_elem(a0, a1); //(z^84 + z^19)u + (z^125 + z^5 + z^3)
	poly64x2_t b0 = {1025, 68719476736};
	poly64x2_t b1 = {2, 137438953473};
	ef_elem b = ef_create_elem(b0, b1); //(z^101 + z^64 + z)u + z^100 + z^10 + 1
	
	//Act
	ef_elem a_inv = ef_inv(a);
	ef_elem b_inv = ef_inv(b);
	ef_elem prod_of_inverses = ef_mull(a_inv, b_inv);
	ef_elem prod = ef_mull(a,b);
	ef_elem inv_of_prod = ef_inv(prod);
	
	//Assert
	uint64_t correct = ef_equal(prod_of_inverses, inv_of_prod);
	assert_true(correct, ctr, "extensionfield: ef_inv_test_prod_of_inverses_is_inverse_of_prod FAILED");
}
	
void ef_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_elem a = ef_rand_elem();
		ef_elem b = ef_rand_elem();
		
		//Act
		ef_elem a_inv = ef_inv(a);
		ef_elem b_inv = ef_inv(b);
		ef_elem prod_of_inverses = ef_mull(a_inv, b_inv);
		ef_elem prod = ef_mull(a,b);
		ef_elem inv_of_prod = ef_inv(prod);
		
		//Assert
		correct = ef_equal(prod_of_inverses, inv_of_prod);
		if(!correct) {
			printf("a: ");
			ef_print_hex_nl(a);
			printf("b: ");
			ef_print_hex_nl(b);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_inv_test_prod_of_inverses_is_inverse_of_prod_rnd FAILED");
}

void ef_inv_test_prod_with_inv_is_one(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {262144, 268435456};
	poly64x2_t a1 = {8589934592, 8192};
	ef_elem a = ef_create_elem(a0, a1); //(z^77 + z^33)u + z^92 + z^18
	poly64x2_t one0 = {1, 0};
	poly64x2_t one1 = {0, 0};
	ef_elem one = ef_create_elem(one0, one1);
	
	//Act
	ef_elem a_inv = ef_inv(a);
	ef_elem a_times_inv = ef_mull(a, a_inv);
	
	//Assert
	uint64_t correct = ef_equal(a_times_inv, one);
	assert_true(correct, ctr, "extensionfield: ef_inv_test_prod_with_inv_is_one FAILED");
}

void ef_inv_test_prod_with_inv_is_one_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ef_elem a = ef_rand_elem();
		poly64x2_t one0 = {1, 0};
		poly64x2_t one1 = {0, 0};
		ef_elem one = ef_create_elem(one0, one1);
		
		//Act
		ef_elem a_inv = ef_inv(a);
		ef_elem a_times_inv = ef_mull(a, a_inv);
		
		//Assert
		correct = ef_equal(a_times_inv, one);
		if(!correct) {
			printf("a: ");
			ef_print_hex_nl(a);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_inv_test_prod_with_inv_is_one_rnd FAILED");
}

void ef_sim_inv_test_crosscheck_inv_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		uint64_t len = 3;
		ef_elem inputs[len];
		ef_elem outputs[len];
		inputs[0] = ef_rand_elem();
		inputs[1] = ef_rand_elem();
		inputs[2] = ef_rand_elem();
		
		//Act
		ef_sim_inv(inputs, outputs, len);
		
		//Assert
		correct = ef_equal(outputs[0], ef_inv(inputs[0])) 
			   && ef_equal(outputs[1], ef_inv(inputs[1])) 
			   && ef_equal(outputs[2], ef_inv(inputs[2]));
		if(!correct) {
			printf("Inputs[0]:\n");
			ef_print_hex_nl(inputs[0]);
			printf("Inputs[1]:\n");
			ef_print_hex_nl(inputs[1]);
			printf("Inputs[2]:\n");
			ef_print_hex_nl(inputs[2]);
			break;
		}
	}
	assert_true(correct, ctr, "extensionfield: ef_sim_inv_test_crosscheck_inv_rnd FAILED");
}

void extensionfield_tests(test_ctr *ctr) {
	ef_create_elem_test_example(ctr);
	
	ef_equal_test_equal(ctr);
	ef_equal_test_notequal(ctr);
	
	ef_add_test_example(ctr);
	ef_add_test_doubling_is_zero(ctr);
	ef_add_test_zero_is_identity(ctr);
	ef_add_test_associative(ctr);
	ef_add_test_commutative(ctr);
	
	ef_mull_test_example(ctr);
	ef_mull_test_associative(ctr);
	ef_mull_test_associative_rnd(ctr);
	ef_mull_test_commutative(ctr);
	ef_mull_test_commutative_rnd(ctr);
	ef_mull_test_one_is_identity(ctr);
	ef_mull_test_zero_is_zero(ctr);
	
	ef_square_test_example(ctr);
	ef_square_test_one_is_one(ctr);
	ef_square_test_zero_is_zero(ctr);
	ef_square_test_crosscheck_pmull(ctr);
	ef_square_test_freshmans_dream(ctr);
	ef_square_test_freshmans_dream_rnd(ctr);
	
	ef_inv_test_example(ctr);
	ef_inv_test_inverse_of_one_is_one(ctr);
	ef_inv_zero_outputs_zero(ctr);
	ef_inv_test_inverse_of_inverse_is_original(ctr);
	ef_inv_test_inverse_of_inverse_is_original_rnd(ctr);
	ef_inv_test_prod_of_inverses_is_inverse_of_prod(ctr);
	ef_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(ctr);
	ef_inv_test_prod_with_inv_is_one(ctr);
	ef_inv_test_prod_with_inv_is_one_rnd(ctr);
	
	ef_sim_inv_test_crosscheck_inv_rnd(ctr);
	
	ef_mull_B_test_example(ctr);
}
