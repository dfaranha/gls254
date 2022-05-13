#include <stdio.h>

#include "ec_scalarmull_tests.h"
#include "../common/ec_scalarmull.h"
#include "../common/ec.h"
#include "../common/extensionfield.h"
#include "../common/basefield.h"
#include "../common/utils.h"

void ec_scalarmull_single_test_example(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj)GEN);
	uint64x2x2_t k = (uint64x2x2_t) {{{1984, 0}, {0, 0}}};

	ef_intrl_elem EX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25)));
	ef_intrl_elem EL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB)));
	ef_intrl_elem EZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_laffine expected = ec_lproj_to_laffine(ec_create_point_lproj(EX, EL, EZ)); // expected = 1984 * GEN

	//Act
	ec_point_laffine actual;
	ec_scalarmull_single_endo_w3_table2D_bulk_ptr(&P, k, &actual);

	//Assert
	uint64_t equal = ec_equal_point_laffine(actual, expected);
	uint64_t on_curve = ec_is_on_curve_laffine(actual);
	printf("%lu\n", on_curve);
	assert_true(equal && on_curve, ctr, "ec: ec_scalar_mull_test_example FAILED");

	uint64_t num_runs = 1;

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		P = ec_rand_point_laffine();

		expected = ec_lproj_to_laffine(ec_scalarmull_single(P, k));
		//Act
		ec_scalarmull_single_endo_w3_table2D_bulk_ptr(&P, k, &actual);

		//Assert
		uint64_t equal = ec_equal_point_laffine(actual, expected);
		uint64_t on_curve = ec_is_on_curve_laffine(actual);

		assert_true(equal && on_curve, ctr, "ec: ec_scalar_mull_test_rand FAILED");
	}
}

void ec_scalarmull_single_test_linearity(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t k = (uint64x2x2_t) {{{7876556743120, 65714569742132121}, {11, 0}}};

	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X51441C4EE272FE55, 0X2DB9775DAEDDE550), bf_create_elem(0X12DD1A65F1D5B480, 0X6CB9034E20AD0EEB)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7AA01DC08C73455A, 0X51F2DF8B2F5FA18C), bf_create_elem(0X6730BC49B9A98F41, 0X57DEBE6DBCE321DE)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X80, 0), bf_create_elem(0X0000000200000000, 0X2000000000)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //(order-1) * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //14329 * GEN

	//Act
	ec_point_laffine add_first = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(ec_add(P, Q)), k);
	ec_point_lproj add_after = ec_add_mixed(ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k), ec_laffine_to_lproj(ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(Q), k)));

	//Assert
	uint64_t correct = ec_equal_point_mixed(add_first, add_after) && ec_is_on_curve_laffine(add_first);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_linearity FAILED");
}

void ec_scalarmull_single_test_negation_order_indifference(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t k = (uint64x2x2_t) {{{76591954, 7695159674259757}, {95124743611111, 56214}}};
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB55AF6853BAA916A, 0X368B1A5434DF7331), bf_create_elem(0X4B3628231E3A83C2, 0XCA5DED000C90027)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //2243156791409331652485 * P

	//Act
	ec_point_laffine neg_first = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(ec_neg(P)), k);
	ec_point_laffine neg_last = ec_neg_laffine(ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k));

	//Assert
	uint64_t correct = ec_equal_point_laffine(neg_first, neg_last);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_negation_order_indifference FAILED");
}

void ec_scalarmull_single_test_order_of_gen_is_order_of_subgroup(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t order_minus_1 = (uint64x2x2_t) SUBGROUP_ORDER;
	order_minus_1.val[0][0] -= 1;

	//Act
	ec_point_laffine gen_inv = ec_scalarmull_single_endo(ec_lproj_to_laffine((ec_point_lproj) GEN), order_minus_1);
	ec_point_laffine gen_inv2 = ec_lproj_to_laffine(ec_scalarmull_single(ec_lproj_to_laffine((ec_point_lproj) GEN), order_minus_1));
	printf("%lu", ec_equal_point_laffine(gen_inv, gen_inv2));
	printf("%lu", ec_is_on_curve(ec_laffine_to_lproj(gen_inv)));
	printf("%lu", ec_is_on_curve(ec_laffine_to_lproj(gen_inv2)));
	ec_print_hex_laffine(gen_inv);
	ec_print_hex_laffine(gen_inv2);

	ec_split_scalar decomp = ec_scalar_decomp(order_minus_1);
	printf("k1 %lu, %lu sign %lu\n", decomp.k1[0], decomp.k1[1], decomp.k1_sign);
	printf("k2 %lu, %lu sign %lu\n", decomp.k2[0], decomp.k2[1], decomp.k2_sign);

	ec_point_lproj gen_plus_gen_inv = ec_add_mixed(gen_inv, (ec_point_lproj) GEN);
	//ec_print_expr(gen_plus_gen_inv);

	//Assert
	uint64_t correct = ec_equal_point_lproj(gen_plus_gen_inv, (ec_point_lproj) INFTY);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_order_of_gen_is_order_of_subgroup FAILED");
}

void ec_scalarmull_single_test_final_k_at_once_same_as_factor_one_at_a_time(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t k1 = (uint64x2x2_t) {{{0XB56483101AD77613, 0X123A67F2366799C}, {0, 0}}};
	uint64x2x2_t k2 = (uint64x2x2_t) {{{0XEEFFEE6752A38174, 0X7544}, {0, 0}}};
	uint64x2x2_t k = (uint64x2x2_t) {{{0XD6D3C77A003A139C, 0XACF7CBD9DC139043}, {0X99A09D4FE2DACA55, 0X85}}}; //k=k1*k2

	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //99921481365893197563 * GEN

	//Act
	ec_point_laffine k1P = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k1);
	ec_point_laffine k1k2P = ec_scalarmull_single_endo_w5_randaccess(k1P, k2);
	ec_point_laffine kP = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k);

	//Assert
	uint64_t correct = ec_equal_point_laffine(kP, k1k2P) && ec_is_on_curve_laffine(k1P) && ec_is_on_curve_laffine(k1k2P);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_final_k_at_once_same_as_factor_one_at_a_time FAILED");
}

void ec_scalarmull_single_test_k_one_is_identity(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t k = (uint64x2x2_t) {{{1, 0}, {0, 0}}};

	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000)));
	ec_point_laffine P = ec_lproj_to_laffine(ec_create_point_lproj(PX, PL, PZ)); //P = 78632917462800214 * GEN

	//Act
	ec_point_laffine actual;
	ec_scalarmull_single_endo_w3_table2D_bulk_ptr(&P, k, &actual);
	//Assert
	uint64_t correct = ec_equal_point_laffine(actual, P) && ec_is_on_curve_laffine(actual);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_k_one_is_identity FAILED");
}

// void ec_scalarmull_double_test_example(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k1 = (uint64x2x2_t) {{{0X4168B72, 0}, {0, 0}}};
// 	uint64x2x2_t k2 = (uint64x2x2_t) {{{0X3F28CB71575AF, 0}, {0, 0}}};
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 133713371337 * GEN
//
// 	ef_elem QX = ef_create_elem(bf_create_elem(0XE20478404AAA1B7F, 0X17A4B025EE434AAD), bf_create_elem(0X3AF8DD6E5CE149C3, 0X4DEC004752AFECF9));
// 	ef_elem QL = ef_create_elem(bf_create_elem(0XA7CFD4415C7C3FE0, 0X142DB23E47B34EF1), bf_create_elem(0X586F58E463F3EEE0, 0X4E8DA1EAC6215EE2));
// 	ef_elem QZ = ef_create_elem(bf_create_elem(0X210, 0), bf_create_elem(0, 0X40000000));
// 	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 949494 * GEN
//
// 	ef_elem EX = ef_create_elem(bf_create_elem(0X54932C0CDE84DD0C, 0X2920B6D886EA380B), bf_create_elem(0XE16091FE1767CA8D, 0X386448A745715798));
// 	ef_elem EL = ef_create_elem(bf_create_elem(0X974AAEB6BCA972D6, 0X23C05EFC539C393B), bf_create_elem(0X80434D5FE9F454A8, 0X1BB6F877D66F45F5));
// 	ef_elem EZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
// 	ec_point_lproj expected = ec_create_point_lproj(EX, EL, EZ);
//
// 	//Act
// 	ec_point_lproj actual = ec_scalarmull_double(P, k1, Q, k2);
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(expected, actual) && ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(actual);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_example FAILED");
// }
//
// void ec_scalarmull_double_test_linearity(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k1 = (uint64x2x2_t) {{{0XABCDEF123, 0XABCDEF123}, {0XABCDEF123, 0XABCDEF123}}};
// 	uint64x2x2_t k2 = (uint64x2x2_t) {{{0XBCADBCADBCAD, 0XFFFFFFFF1}, {0X77CC77, 0XAAA}}};
//
// 	ef_elem P1X = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
// 	ef_elem P1L = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
// 	ef_elem P1Z = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
// 	ec_point_lproj P1 = ec_create_point_lproj(P1X, P1L, P1Z); //P1 = 133713371337 * GEN
//
// 	ef_elem P2X = ef_create_elem(bf_create_elem(0XE232234BCDA1F7A9, 0X72E20184752C18C3), bf_create_elem(0XA249666BD031EF41,0X78BFF5D2CAFC6FC2));
// 	ef_elem P2L = ef_create_elem(bf_create_elem(0X1064B65732340CD9, 0X63C8BCD9FBA1773F), bf_create_elem(0X128743886A3EC579,0X6E8E5EC1D3091E7F));
// 	ef_elem P2Z = ef_create_elem(bf_create_elem(0X212FB30BAC886248, 0X797A835808F9BF57), bf_create_elem(0X8CA3223B872913CE,0X4B0C008D15FE6EB1));
// 	ec_point_lproj P2 = ec_create_point_lproj(P2X, P2L, P2Z); //P2 = 78632917462800214 * 2 * GEN
//
// 	ef_elem P3X = ef_create_elem(bf_create_elem(0XE20478404AAA1B7F, 0X17A4B025EE434AAD), bf_create_elem(0X3AF8DD6E5CE149C3, 0X4DEC004752AFECF9));
// 	ef_elem P3L = ef_create_elem(bf_create_elem(0XA7CFD4415C7C3FE0, 0X142DB23E47B34EF1), bf_create_elem(0X586F58E463F3EEE0, 0X4E8DA1EAC6215EE2));
// 	ef_elem P3Z = ef_create_elem(bf_create_elem(0X210, 0), bf_create_elem(0, 0X40000000));
// 	ec_point_lproj P3 = ec_create_point_lproj(P3X, P3L, P3Z); //P3 = 949494 * GEN
//
// 	ef_elem P4X = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
// 	ef_elem P4L = ef_create_elem(bf_create_elem(0XB55AF6853BAA91FA, 0X368B1A5434DF7331), bf_create_elem(0X4B36A8231E3A83C2, 0XCA5DED000C90427));
// 	ef_elem P4Z = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
// 	ec_point_lproj P4 = ec_create_point_lproj(P4X, P4L, P4Z); //P4 = -2243156791409331652485 * GEN
//
// 	//Act
// 	ec_point_lproj result_double = ec_scalarmull_double(ec_add(P1, P2), k1, ec_add(P3, P4), k2);
// 	ec_point_lproj result_single = ec_add(ec_add(ec_add(ec_scalarmull_single_endo_w5_randaccess(P1, k1), ec_scalarmull_single_endo_w5_randaccess(P2, k1)), ec_scalarmull_single_endo_w5_randaccess(P3, k2)), ec_scalarmull_single_endo_w5_randaccess(P4, k2));
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(result_double, result_single) && ec_is_on_curve(result_double);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_linearity FAILED");
// }
//
// void ec_scalarmull_double_test_commutative(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k1 = (uint64x2x2_t) {{{0XAAA, 0XBBB}, {0XCCC, 0XDDD}}};
// 	uint64x2x2_t k2 = (uint64x2x2_t) {{{0XEEE, 0XFFF}, {0X111, 0X222}}};
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 133713371337 * GEN
//
// 	ef_elem QX = ef_create_elem(bf_create_elem(0XE232234BCDA1F7A9, 0X72E20184752C18C3), bf_create_elem(0XA249666BD031EF41,0X78BFF5D2CAFC6FC2));
// 	ef_elem QL = ef_create_elem(bf_create_elem(0X1064B65732340CD9, 0X63C8BCD9FBA1773F), bf_create_elem(0X128743886A3EC579,0X6E8E5EC1D3091E7F));
// 	ef_elem QZ = ef_create_elem(bf_create_elem(0X212FB30BAC886248, 0X797A835808F9BF57), bf_create_elem(0X8CA3223B872913CE,0X4B0C008D15FE6EB1));
// 	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 78632917462800214 * 2 * GEN
//
// 	//Act
// 	ec_point_lproj k1P_plus_k2Q = ec_scalarmull_double(P, k1, Q, k2);
// 	ec_point_lproj k2Q_plus_k1P = ec_scalarmull_double(Q, k2, P, k1);
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(k1P_plus_k2Q, k2Q_plus_k1P) && ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(k1P_plus_k2Q);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_commutative FAILED");
// }
//
// void ec_scalarmull_double_test_negation_order_indifference(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k1 = (uint64x2x2_t) {{{0XA7A7A, 0XB8B8B}, {0XC9C9C, 0XD0D0D}}};
// 	uint64x2x2_t k2 = (uint64x2x2_t) {{{0XE7E7E, 0XF8F8F}, {0X19191, 0X20202}}};
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //78632917462800214 * GEN
//
// 	ef_elem QX = ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25));
// 	ef_elem QL = ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB));
// 	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
// 	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , z = u + z
//
// 	//Act
// 	ec_point_lproj neg_first = ec_scalarmull_double(ec_neg(P), k1, ec_neg(Q), k2);
// 	ec_point_lproj neg_after = ec_neg(ec_scalarmull_double(P, k1, Q, k2));
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(neg_first, neg_after) && ec_is_on_curve(neg_first);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_negation_order_indifference FAILED");
// }
//
// void ec_scalarmull_double_test_k_ones_is_add(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k = (uint64x2x2_t) {{{1, 0}, {0, 0}}};
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 99921481365893197563 * GEN
//
// 	ef_elem QX = ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78));
// 	ef_elem QL = ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF));
// 	ef_elem QZ = ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796));
// 	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 14329 * GEN
//
// 	//Act
// 	ec_point_lproj scalarmull_result = ec_scalarmull_double(P, k, Q, k);
// 	ec_point_lproj add_result = ec_add(P, Q);
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(scalarmull_result, add_result);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_k_ones_is_add FAILED");
// }
//
// void ec_scalarmull_double_test_crosscheck_scalarmull_single(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k1 = (uint64x2x2_t) {{{12345, 6789}, {10111213, 141516}}};
// 	uint64x2x2_t k2 = (uint64x2x2_t) {{{24690, 13578}, {20222426, 283032}}}; //k2 = 2 * k1
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0XB55AF6853BAA91FA, 0X368B1A5434DF7331), bf_create_elem(0X4B36A8231E3A83C2, 0XCA5DED000C90427));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = -2243156791409331652485 * GEN
//
// 	//Act
// 	ec_point_lproj double_result = ec_scalarmull_double(P, k1, P, k1);
// 	ec_point_lproj single_result = ec_scalarmull_single_endo_w5_randaccess(P, k2);
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(double_result, single_result);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_crosscheck_scalarmull_single FAILED");
// }
//
// void ec_scalarmull_double_test_point_and_neg_cancel(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k = (uint64x2x2_t) {{{999213654, 2084685670}, {32123142, 6462462412477}}};
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 133713371337 * GEN
//
// 	//Act
// 	ec_point_lproj result = ec_scalarmull_double(P, k, ec_neg(P), k);
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj((ec_point_lproj) INFTY, result);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_point_and_neg_cancel FAILED");
// }
//
// void ec_scalarmull_double_test_point_and_neg_interfere(test_ctr *ctr) {
// 	//Arrange
// 	uint64x2x2_t k1 = (uint64x2x2_t) {{{1, 0}, {0, 0}}};
// 	uint64x2x2_t k2 = (uint64x2x2_t) {{{2, 0}, {0, 0}}};
//
// 	ef_elem PX = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
// 	ef_elem PL = ef_create_elem(bf_create_elem(0XB55AF6853BAA916A, 0X368B1A5434DF7331), bf_create_elem(0X4B3628231E3A83C2, 0XCA5DED000C90027));
// 	ef_elem PZ = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
// 	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //2243156791409331652485 * P
//
// 	//Act
// 	ec_point_lproj P_neg = ec_neg(P);
// 	ec_point_lproj result = ec_scalarmull_double(P, k1, P_neg, k2);
//
// 	//Assert
// 	uint64_t correct = ec_equal_point_lproj(P_neg, result);
// 	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_point_and_neg_interfere FAILED");
// }

void ec_scalarmull_single_endo_test_example(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t k = (uint64x2x2_t) {{{12345, 0}, {0,0}}};

	ec_point_lproj P_lproj = (ec_point_lproj) GEN;
	ec_point_laffine P_laffine = ec_lproj_to_laffine(P_lproj);

	ef_intrl_elem EX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem EL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ec_point_laffine E = ec_create_point_laffine(EX, EL); //E = 12345 * GEN

	//Act
	ec_point_laffine R = ec_scalarmull_single_endo(P_laffine, k);

	//Assert
	uint64_t correct = ec_equal_point_laffine(R, E);
	assert_true(correct, ctr, "ec_scalarmull_single_endo_test_example FAILED");
}

void ec_scalarmull_single_endo_test_crosscheck_rnd(test_ctr *ctr) {
	//Arrange
	uint64x2x2_t k = ec_rand_scalar();
	ec_point_laffine P = ec_rand_point_laffine();

	//Act
	ec_point_laffine kP = ec_scalarmull_single_endo_w5_randaccess(P, k);
	ec_point_laffine crosscheck = ec_scalarmull_single_endo(P, k);

	//Assert
	uint64_t equal = ec_equal_point_laffine(kP, crosscheck);
	if(!equal) {
		printf("k: %lu, %lu, %lu, %lu\n", k.val[0][0], k.val[0][1], k.val[1][0], k.val[1][1]);
		ec_print_hex_laffine(P);
	}
	assert_true(equal, ctr, "ec_scalarmull_single_endo_test_crosscheck_rnd FAILED");
}

void ec_scalarmull_single_rand_test_w6(test_ctr *ctr) {
	uint64x2x2_t k = ec_rand_scalar();
	ec_point_laffine P = ec_rand_point_laffine();

	ec_point_lproj expected = ec_scalarmull_single(P, k);
	//Act
	ec_point_laffine actual = ec_scalarmull_single_endo_w6_randaccess(P, k);

	//Assert
	uint64_t equal = ec_equal_point_mixed(actual, expected);
	uint64_t on_curve = ec_is_on_curve_laffine(actual);

	assert_true(equal && on_curve, ctr, "ec: ec_scalarmull_single_rand_test_w6 FAILED");
}

void ec_scalarmull_single_rand_test_w4(test_ctr *ctr) {
	uint64x2x2_t k = ec_rand_scalar();
	ec_point_laffine P = ec_rand_point_laffine();

	ec_point_lproj expected = ec_scalarmull_single(P, k);
	//Act
	ec_point_laffine actual = ec_scalarmull_single_endo_w4_randaccess(P, k);

	//Assert
	uint64_t equal = ec_equal_point_mixed(actual, expected);
	uint64_t on_curve = ec_is_on_curve_laffine(actual);

	assert_true(equal && on_curve, ctr, "ec: ec_scalarmull_single_rand_test_w4 FAILED");
}

void ec_scalarmull_single_rand_test_w3(test_ctr *ctr) {
	uint64x2x2_t k = ec_rand_scalar();
	ec_point_laffine P = ec_rand_point_laffine();

	ec_point_lproj expected = ec_scalarmull_single(P, k);
	//Act
	ec_point_laffine actual = ec_scalarmull_single_endo_w3_randaccess(P, k);

	//Assert
	uint64_t equal = ec_equal_point_mixed(actual, expected);
	uint64_t on_curve = ec_is_on_curve_laffine(actual);

	assert_true(equal && on_curve, ctr, "ec: ec_scalarmull_single_rand_test_w3 FAILED");
}

void ec_scalarmull_test_precomputation(test_ctr *ctr) {
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //99921481365893197563 * GEN

	ec_point_laffine table[16];
	precompute_w5(ec_lproj_to_laffine(P), table);

	uint64_t equal1 = ec_equal_point_laffine(table[1], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{1, 0}, {0, 0}}}));
	uint64_t equal3 = ec_equal_point_laffine(table[3], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{3, 0}, {0, 0}}}));
	uint64_t equal5 = ec_equal_point_laffine(table[5], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{5, 0}, {0, 0}}}));
	uint64_t equal7 = ec_equal_point_laffine(table[7], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{7, 0}, {0, 0}}}));
	uint64_t equal9 = ec_equal_point_laffine(table[9], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{9, 0}, {0, 0}}}));
	uint64_t equal11 = ec_equal_point_laffine(table[11], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{11, 0}, {0, 0}}}));
	uint64_t equal13 = ec_equal_point_laffine(table[13], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{13, 0}, {0, 0}}}));
	uint64_t equal15 = ec_equal_point_laffine(table[15], ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), (uint64x2x2_t) {{{15, 0}, {0, 0}}}));

	assert_true(equal1 && equal3 && equal5 && equal7 && equal9 && equal11 && equal13 && equal15 , ctr, "ec_scalarmull_test_precomputation FAILED");
}

void ec_scalarmull_test_precomputation_w4_table2D_ptr(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj) GEN);
	ec_point_laffine table[16];
	
	//Act
	ec_precompute_w4_table2D_ptr(&P, table);
	
	//Assert
	uint64_t correct = 1;
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {
			uint64x2x2_t k1 = (uint64x2x2_t) {{{2*i+1, 0}, {0, 0}}};
			uint64x2x2_t k2 = (uint64x2x2_t) {{{2*j+1, 0}, {0, 0}}};
			ec_point_laffine k1P = ec_scalarmull_single_endo_w5_randaccess(P, k1);
			ec_point_laffine k2EndoP = ec_scalarmull_single_endo_w5_randaccess(ec_endo_laffine(P), k2);
			ec_point_laffine expected = ec_lproj_to_laffine(ec_add_laffine_unchecked(k1P, k2EndoP));
			uint64_t equal = ec_equal_point_laffine(expected, table[4*i+j]);
			correct &= equal; 
			if(!equal) {
				printf("i: %d j: %d\n", i, j);
			}
		}
	}
	assert_true(correct, ctr, "ec_scalarmull_test_precomputation_w4_table2D FAILED");
}

void ec_scalarmull_test_precomputation_w3_table2D_ptr(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj) GEN);
	ec_point_laffine table[4];
	
	//Act
	ec_precompute_w3_table2D_ptr(&P, table);
	
	//Assert
	uint64_t correct = 1;
	for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 2; i++) {
			uint64x2x2_t k1 = (uint64x2x2_t) {{{2*i+1, 0}, {0, 0}}};
			uint64x2x2_t k2 = (uint64x2x2_t) {{{2*j+1, 0}, {0, 0}}};
			ec_point_laffine k1P = ec_scalarmull_single_endo_w5_randaccess(P, k1);
			ec_point_laffine k2EndoP = ec_scalarmull_single_endo_w5_randaccess(ec_endo_laffine(P), k2);
			ec_point_laffine expected = ec_lproj_to_laffine(ec_add_laffine_unchecked(k1P, k2EndoP));
			uint64_t equal = ec_equal_point_laffine(expected, table[2*i+j]);
			correct &= equal; 
			if(!equal) {
				printf("i: %d j: %d\n", i, j);
			}
		}
	}
	assert_true(correct, ctr, "ec_scalarmull_test_precomputation_w3_table2D FAILED");
}

void ec_scalarmull_test_lookup_from_w4_table2D_ptr(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj) GEN);
	ec_point_laffine table[16];
	ec_precompute_w4_table2D_ptr(&P, table);
	ec_split_scalar decomp;
	decomp.k1_sign = 0;
	decomp.k2_sign = 0;
	signed char rec_k1[1];
	signed char rec_k2[1];
	signed char k1 = 3;
	signed char k2 = 5;
	
	ec_point_laffine iP = ec_scalarmull_single_endo_w5_randaccess(P, (uint64x2x2_t) {{{k1, 0}, {0, 0}}});
	ec_point_laffine jEndoP = ec_scalarmull_single_endo_w5_randaccess(ec_endo_laffine(P), (uint64x2x2_t) {{{k2, 0}, {0, 0}}});
	ec_point_laffine iPneg = ec_neg_laffine(iP);
	ec_point_laffine jEndoPneg = ec_neg_laffine(jEndoP);

	rec_k1[0] = k1;
	rec_k2[0] = k2;
	ec_point_lproj tmp;
	ec_add_laffine_unchecked_ptr(&iP, &jEndoP, &tmp);
	ec_point_laffine expected = ec_lproj_to_laffine(tmp);
	ec_point_laffine actual;
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	uint64_t equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w4_table2D_ptr +i +j FAILED");

	rec_k1[0] = k1;
	rec_k2[0] = -k2;
	ec_add_laffine_unchecked_ptr(&iP, &jEndoPneg, &tmp);
	expected = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w4_table2D_ptr +i -j FAILED");

	rec_k1[0] = -k1;
	rec_k2[0] = k2;
	ec_add_laffine_unchecked_ptr(&iPneg, &jEndoP, &tmp);
	expected = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w4_table2D_ptr -i +j FAILED");

	rec_k1[0] = -k1;
	rec_k2[0] = -k2;
	ec_add_laffine_unchecked_ptr(&iPneg, &jEndoPneg, &tmp);
	expected = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w4_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w4_table2D_ptr -i -j FAILED");
}

void ec_scalarmull_test_lookup_from_w3_table2D_ptr(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj) GEN);
	ec_point_laffine table[4];
	ec_precompute_w3_table2D_ptr(&P, table);
	ec_split_scalar decomp;
	decomp.k1_sign = 0;
	decomp.k2_sign = 0;
	signed char rec_k1[1];
	signed char rec_k2[1];
	signed char k1 = 3;
	signed char k2 = 1;
	
	ec_point_laffine iP = ec_scalarmull_single_endo_w5_randaccess(P, (uint64x2x2_t) {{{k1, 0}, {0, 0}}});
	ec_point_laffine jEndoP = ec_endo_laffine(P);
	ec_point_laffine iPneg = ec_neg_laffine(iP);
	ec_point_laffine jEndoPneg = ec_neg_laffine(jEndoP);

	rec_k1[0] = k1;
	rec_k2[0] = k2;
	ec_point_lproj tmp;
	ec_add_laffine_unchecked_ptr(&iP, &jEndoP, &tmp);
	ec_point_laffine expected = ec_lproj_to_laffine(tmp);
	ec_point_laffine actual;
	ec_lookup_from_w3_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	printf("%lu", ec_is_on_curve_laffine(actual));
	uint64_t equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w3_table2D_ptr +i +j FAILED");

	rec_k1[0] = k1;
	rec_k2[0] = -k2;
	ec_add_laffine_unchecked_ptr(&iP, &jEndoPneg, &tmp);
	expected = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w3_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w3_table2D_ptr +i -j FAILED");

	rec_k1[0] = -k1;
	rec_k2[0] = k2;
	ec_add_laffine_unchecked_ptr(&iPneg, &jEndoP, &tmp);
	expected = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w3_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w3_table2D_ptr -i +j FAILED");

	rec_k1[0] = -k1;
	rec_k2[0] = -k2;
	ec_add_laffine_unchecked_ptr(&iPneg, &jEndoPneg, &tmp);
	expected = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w3_table2D_ptr(&decomp, rec_k1, rec_k2, table, 0, &actual);
	equal = ec_equal_point_laffine(expected, actual);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w3_table2D_ptr -i -j FAILED");
}

void ec_scalarmull_test_lookup_from_w3_table2D_bulk_ptr(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj) GEN);
	ec_point_laffine table[4];
	ec_precompute_w3_table2D_ptr(&P, table);
	ec_split_scalar decomp;
	decomp.k1_sign = 0;
	decomp.k2_sign = 0;
	signed char rec_k1[2];
	signed char rec_k2[2];
	signed char k1 = 1;
	signed char k2 = 3;

	ec_point_laffine iP = P;
	ec_point_laffine jEndoP = ec_scalarmull_single_endo_w5_randaccess(ec_endo_laffine(P), (uint64x2x2_t) {{{k2, 0}, {0, 0}}});
	ec_point_laffine iPneg = ec_neg_laffine(iP);
	ec_point_laffine jEndoPneg = ec_neg_laffine(jEndoP);

	rec_k1[0] = k1;
	rec_k2[0] = k2;
	rec_k1[1] = -k1;
	rec_k2[1] = -k2;
	ec_point_lproj tmp;
	ec_add_laffine_unchecked_ptr(&iP, &jEndoP, &tmp);
	ec_point_laffine expected1 = ec_lproj_to_laffine(tmp);
	ec_add_laffine_unchecked_ptr(&iPneg, &jEndoPneg, &tmp);
	ec_point_laffine expected2 = ec_lproj_to_laffine(tmp);
	ec_point_laffine actual1, actual2;
	ec_lookup_from_w3_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, 0, 1, &actual1, &actual2);
	uint64_t equal = ec_equal_point_laffine(expected1, actual1) && ec_equal_point_laffine(expected2, actual2);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w3_table2D_bulk_ptr same sign FAILED");

	rec_k1[0] = k1;
	rec_k2[0] = -k2;
	rec_k1[1] = -k1;
	rec_k2[1] = k2;
	ec_add_laffine_unchecked_ptr(&iP, &jEndoPneg, &tmp);
	expected1 = ec_lproj_to_laffine(tmp);
	ec_add_laffine_unchecked_ptr(&iPneg, &jEndoP, &tmp);
	expected2 = ec_lproj_to_laffine(tmp);
	ec_lookup_from_w3_table2D_bulk_ptr(&decomp, rec_k1, rec_k2, table, 0, 1, &actual1, &actual2);
	equal = ec_equal_point_laffine(expected1, actual1) && ec_equal_point_laffine(expected2, actual2);
	assert_true(equal, ctr, "ec_scalarmull_test_lookup_from_w3_table2D_ptr neq sign FAILED");
}

void ec_scalarmull_tests(test_ctr *ctr) {
	ec_scalarmull_single_test_example(ctr);
	/*ec_scalarmull_single_test_linearity(ctr);
	ec_scalarmull_single_test_negation_order_indifference(ctr);
	ec_scalarmull_single_test_order_of_gen_is_order_of_subgroup(ctr);
	ec_scalarmull_single_test_final_k_at_once_same_as_factor_one_at_a_time(ctr);*/
	ec_scalarmull_single_test_k_one_is_identity(ctr);

	// ec_scalarmull_double_test_example(ctr);
	// ec_scalarmull_double_test_linearity(ctr);
	// ec_scalarmull_double_test_commutative(ctr);
	// ec_scalarmull_double_test_negation_order_indifference(ctr);
	// ec_scalarmull_double_test_k_ones_is_add(ctr);
	// ec_scalarmull_double_test_crosscheck_scalarmull_single(ctr);
	// ec_scalarmull_double_test_point_and_neg_cancel(ctr);
	// ec_scalarmull_double_test_point_and_neg_interfere(ctr);

	/*ec_scalarmull_single_endo_test_example(ctr);
	ec_scalarmull_single_endo_test_crosscheck_rnd(ctr);

	ec_scalarmull_single_rand_test_w6(ctr);
	ec_scalarmull_single_rand_test_w4(ctr);
	ec_scalarmull_single_rand_test_w3(ctr);*/

	// ec_scalarmull_test_precomputation(ctr);
	ec_scalarmull_test_precomputation_w4_table2D_ptr(ctr);
	ec_scalarmull_test_precomputation_w3_table2D_ptr(ctr);
	ec_scalarmull_test_lookup_from_w4_table2D_ptr(ctr);
	ec_scalarmull_test_lookup_from_w3_table2D_ptr(ctr);
	ec_scalarmull_test_lookup_from_w3_table2D_bulk_ptr(ctr);
}
