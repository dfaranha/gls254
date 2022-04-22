#include "util_tests.h"
#include "../common/utils.h"

void compare_doubles_test_equal(test_ctr *ctr) {
	//Arrange
	double a = 9.3;
	double b = 9.3;
	double errmargin = 0.0000001;
	
	//Act
	uint64_t are_equal = compare_doubles(a, b, errmargin) == 0;
	
	//Assert
	assert_true(are_equal, ctr, "utils: compare_doubles_test_equal FAILED");
}

void compare_doubles_test_approxequal(test_ctr *ctr) {
	//Arrange
	double a = -9.3;
	double b = -9.30000001;
	double errmargin = 0.0000001;
	
	//Act
	uint64_t are_equal = compare_doubles(a, b, errmargin) == 0;
	
	//Assert
	assert_true(are_equal, ctr, "utils: compare_doubles_test_approxequal FAILED");
}

void compare_doubles_test_notequal(test_ctr *ctr) {
	//Arrange
	double a = 9.3;
	double b = 9.3000002;
	double errmargin = 0.0000001;
	
	//Act & Assert
	uint64_t are_equal = compare_doubles(a, b, errmargin) == 0;
	assert_false(are_equal, ctr, "utils: compare_doubles_test_notequal FAILED");
}

void compare_doubles_test_lessthan(test_ctr *ctr) {
	//Arrange
	double a = 9.3;
	double b = 9.3000002;
	double errmargin = 0.0000001;
	
	//Act
	uint64_t b_greater = compare_doubles(a, b, errmargin) == -1;
	
	//Assert
	assert_true(b_greater, ctr, "utils: compare_doubles_test_lessthan FAILED");
}

void compare_doubles_test_greaterthan(test_ctr *ctr) {
	//Arrange
	double a = -13.111111;
	double b = -13.111112;
	double errmargin = 0.0000001;
	
	//Act
	uint64_t a_greater = compare_doubles(a, b, errmargin) == 1;
	
	//Assert
	assert_true(a_greater, ctr, "utils: compare_doubles_test_greaterthan FAILED");
}

void equal_poly64x2_test_equal(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {123, 45678910};
	poly64x2_t b = {123, 45678910};
	
	//Act
	uint64_t equal = equal_poly64x2(a, b);
	
	//Assert
	assert_true(equal, ctr, "utils: equal_poly64x2_test_equal FAILED");
}

void equal_poly64x2_test_notequal(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {123, 45678910};
	poly64x2_t b = {123, 45688910};
	
	//Act
	uint64_t equal = equal_poly64x2(a, b);
	
	//Assert
	assert_false(equal, ctr, "utils: equal_poly64x2_test_notequal FAILED");
}

void equal_poly64x2x2_test_equal(test_ctr *ctr) {
	//Arrange
	poly64x2x2_t a;
	poly64x2_t a0 = {8594687, 97466431};
	poly64x2_t a1 = {93651237, 10389541};
	a.val[0] = a0;
	a.val[1] = a1;
	poly64x2x2_t b;
	poly64x2_t b0 = {8594687, 97466431};
	poly64x2_t b1 = {93651237, 10389541};
	b.val[0] = b0;
	b.val[1] = b1;
	
	//Act
	uint64_t equal = equal_poly64x2x2(a, b);
	
	//Assert
	assert_true(equal, ctr, "utils: equal_poly64x2x2_test_equal FAILED");
}

void equal_poly64x2x2_test_notequal(test_ctr *ctr) {
	//Arrange
	poly64x2x2_t a;
	poly64x2_t a0 = {8594687, 97466431};
	poly64x2_t a1 = {93651237, 10389541};
	a.val[0] = a0;
	a.val[1] = a1;
	poly64x2x2_t b;
	poly64x2_t b0 = {8594686, 97466431};
	poly64x2_t b1 = {93651237, 10389541};
	b.val[0] = b0;
	b.val[1] = b1;
	
	//Act
	uint64_t equal = equal_poly64x2x2(a, b);
	
	//Assert
	assert_false(equal, ctr, "utils: equal_poly64x2x2_test_notequal FAILED");
}

void concat_bf_poly_test(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {76194397641, 9487521};
	poly64x2_t a1 = {333666999, 71421283542};
	
	//Act
	poly64x2x2_t a = concat_bf_poly(a0, a1);
	
	//Assert
	uint64_t equal = equal_poly64x2(a.val[0], a0) && equal_poly64x2(a.val[1], a1);
	assert_true(equal, ctr, "utils: concat_bf_poly_test FAILED");
}

void average_test(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = { 15, 2, 9};
	double expected = (15.0+2.0+9.0)/3.0;
	double errmargin = 0.0000001;
	
	//Act
	double actual = average(times, 3);
	
	//Assert
	uint64_t in_err_interval = compare_doubles(expected, actual, errmargin) == 0;
	assert_true(in_err_interval, ctr, "utils: average_test FAILED");
}

void median_test_even(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {2, 15, 17, 18, 23, 36};
	double expected = 17.5;
	double errmargin = 0.0000001;
	
	//Act
	double actual = median(times, 6);
	
	//Assert
	uint64_t in_err_interval = compare_doubles(expected, actual, errmargin) == 0;
	assert_true(in_err_interval, ctr, "utils: median_test_even FAILED");
}

void median_test_odd(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {2, 15, 17, 18, 20, 23, 36};
	double expected = 18;
	double errmargin = 0.0000001;
	
	//Act
	double actual = median(times, 7);
	
	//Assert
	uint64_t in_err_interval = compare_doubles(expected, actual, errmargin) == 0;
	assert_true(in_err_interval, ctr, "utils: median_test_odd FAILED");
}

void util_tests(test_ctr *ctr) {
	compare_doubles_test_equal(ctr);
	compare_doubles_test_approxequal(ctr);
	compare_doubles_test_notequal(ctr);
	compare_doubles_test_lessthan(ctr);
	compare_doubles_test_greaterthan(ctr);
	equal_poly64x2_test_equal(ctr);
	equal_poly64x2_test_notequal(ctr);
	equal_poly64x2x2_test_equal(ctr);
	equal_poly64x2x2_test_notequal(ctr);
	concat_bf_poly_test(ctr);
	average_test(ctr);
	median_test_even(ctr);
	median_test_odd(ctr);
}
