#include <stdio.h>

#include "ec_tests.h"
#include "../common/ec.h"
#include "../common/utils.h"

void ec_create_point_lproj_test_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem x = ef_intrl_interleave(ef_create_elem(bf_create_elem(16785442, 721476), bf_create_elem(78099554548664, 6547959942615)));
	ef_intrl_elem l = ef_intrl_interleave(ef_create_elem(bf_create_elem(972573, 353523623), bf_create_elem(23523988509283, 2435335)));
	ef_intrl_elem z = ef_intrl_interleave(ef_create_elem(bf_create_elem(54689374, 4584), bf_create_elem(34548734, 5648574685)));

	//Act
	ec_point_lproj P = ec_create_point_lproj(x, l, z);

	//Assert
	uint64_t correct = ef_intrl_equal(P.x, x) && ef_intrl_equal(P.l, l) && ef_intrl_equal(P.z, z);
	assert_true(correct, ctr, "ec: ec_create_point_lproj_test_example FAILED");
}

void ec_equal_point_lproj_test_equivalent_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(8580737671568810207U, 2871255672283442416U), bf_create_elem(14038621212797386049U, 7102795227656941388U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(8047436918354421319U, 1946646612861660792U), bf_create_elem(3809596628914525439U, 6366755232823307697U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0,0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 1984 * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , Q.z = u + z

	//Act & Assert
	uint64_t equal = ec_equal_point_lproj(P, Q);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q);
	assert_true(equal && on_curve, ctr, "ec: ec_equal_point_lproj_test_equivalent_example FAILED");
}

void ec_equal_point_lproj_test_notequivalent_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , Q.z = u + z

	//Act & Assert
	uint64_t equal = ec_equal_point_lproj(P, Q);
	assert_false(equal, ctr, "ec: ec_equal_point_lproj_test_notequivalent_example FAILED");
}

void ec_is_on_curve_test_generator_on_curve(test_ctr *ctr) {
	//Arrange, Act & Assert
	uint64_t on_curve = ec_is_on_curve((ec_point_lproj) GEN);
	assert_true(on_curve, ctr, "ec: test_generator_on_curve FAILED");
}

void ec_is_on_curve_test_point_not_on_curve(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437212U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ);

	//Act & Assert
	uint64_t on_curve = ec_is_on_curve((ec_point_lproj) P);
	assert_false(on_curve, ctr, "ec: test_generator_on_curve FAILED");
}


void ec_rand_point_lproj_test_on_curve(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange & Act
		ec_point_lproj P = ec_rand_point_lproj();

		//Assert
		correct = ec_is_on_curve(P);
		if(!correct) {
			printf("p: ");
			ec_print_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_rand_point_lproj_test_on_curve FAILED");
}

void ec_rand_point_laffine_test_on_curve(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange & Act
		ec_point_laffine P = ec_rand_point_laffine();
		ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};
		ec_point_lproj L;

		L.x = P.x;
		L.l = P.l;
		L.z = one;


		//Assert
		correct = ec_is_on_curve(L);
		if(!correct) {
			printf("p: ");
			ec_print_hex(L);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_rand_point_laffine_test_on_curve FAILED");
}

void ec_add_test_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , z = u + z

	ef_intrl_elem eX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78)));
	ef_intrl_elem eL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF)));
	ef_intrl_elem eZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796)));
	ec_point_lproj expected = ec_create_point_lproj(eX, eL, eZ);

	//Act
	ec_point_lproj actual = ec_add(P, Q);

	//Assert
	uint64_t equal = ec_equal_point_lproj(expected, actual);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(actual);
	assert_true(equal && on_curve, ctr, "ec: ec_add_test_example FAILED");
}

void ec_add_test_is_on_curve_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj Q = ec_rand_point_lproj();
		
		//Act
		ec_point_lproj sum = ec_add(P, Q);

		//Assert
		correct = ec_is_on_curve(sum);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			printf("Q: ");
			ec_print_hex(Q);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_add_test_is_on_curve_rnd FAILED");
}

void ec_add_test_associative(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(5130258657669722291U, 2683433950625433362U), bf_create_elem(15861652668403718055U, 2280238350963310U)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(13076644468273807311U, 8504358646598259325U), bf_create_elem(16381575749271831681U, 1938714279827322046U)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 674559848943297 * GEN

	ef_intrl_elem RX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25)));
	ef_intrl_elem RL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB)));
	ef_intrl_elem RZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_lproj R = ec_create_point_lproj(RX, RL, RZ); // Q = 1984 * GEN , z = u + z

	//Act
	ec_point_lproj P_plus_Q = ec_add(P, Q);
	ec_point_lproj P_plus_Q_first = ec_add(P_plus_Q, R);
	ec_point_lproj Q_plus_R = ec_add(Q, R);
	ec_point_lproj Q_plus_R_first = ec_add(P, Q_plus_R);

	//Assert
	uint64_t equal = ec_equal_point_lproj(P_plus_Q_first, Q_plus_R_first);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(R) && ec_is_on_curve(P_plus_Q) && ec_is_on_curve(P_plus_Q_first) && ec_is_on_curve(Q_plus_R);
	assert_true(equal && on_curve, ctr, "ec: ec_add_test_associative FAILED");
}

void ec_add_test_associative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj Q = ec_rand_point_lproj();
		ec_point_lproj R = ec_rand_point_lproj();

		//Act
		ec_point_lproj P_plus_Q = ec_add(P, Q);
		ec_point_lproj P_plus_Q_first = ec_add(P_plus_Q, R);
		ec_point_lproj Q_plus_R = ec_add(Q, R);
		ec_point_lproj Q_plus_R_first = ec_add(P, Q_plus_R);

		//Assert
		correct = ec_equal_point_lproj(P_plus_Q_first, Q_plus_R_first);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			printf("Q: ");
			ec_print_hex(Q);
			printf("R: ");
			ec_print_hex(R);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_add_test_associative_rnd FAILED");
}

void ec_add_test_commutative(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , z = u + z

	//Act
	ec_point_lproj P_plus_Q = ec_add(P, Q);
	ec_point_lproj Q_plus_P = ec_add(Q, P);

	//Assert
	uint64_t equal = ec_equal_point_lproj(P_plus_Q, Q_plus_P);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(P_plus_Q);
	assert_true(equal && on_curve, ctr, "ec: ec_add_test_commutative FAILED");
}

void ec_add_test_commutative_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj Q = ec_rand_point_lproj();

		//Act
		ec_point_lproj P_plus_Q = ec_add(P, Q);
		ec_point_lproj Q_plus_P = ec_add(Q, P);

		//Assert
		correct = ec_equal_point_lproj(P_plus_Q, Q_plus_P);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			printf("Q: ");
			ec_print_hex(Q);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_add_test_commutative_rnd FAILED");
}

//Might want to remove these three for future runtime optimizations?
void ec_add_test_point_plus_infty_is_point(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); // P = 1984 * GEN , z = u + z

	//Act
	ec_point_lproj result = ec_add(P, (ec_point_lproj) INFTY);

	//Assert
	uint64_t correct = ec_equal_point_lproj(P, result);
	assert_true(correct, ctr, "ec: ec_add_test_point_plus_infty_is_point FAILED");
}

void ec_add_test_infty_plus_point_is_point(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //12345 * GEN + 1984 * GEN

	//Act
	ec_point_lproj result = ec_add((ec_point_lproj) INFTY, P);

	//Assert
	uint64_t correct = ec_equal_point_lproj(P, result);
	assert_true(correct, ctr, "ec: ec_add_test_point_plus_infty_is_point FAILED");
}

void ec_add_test_crosscheck_double(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //12345 * GEN + 1984 * GEN

	//Act
	ec_point_lproj sum = ec_add(P, P);
	ec_point_lproj doubling = ec_double(P);

	//Assert
	uint64_t correct = ec_equal_point_lproj(sum, doubling) && ec_is_on_curve(P) && ec_is_on_curve(sum);
	assert_true(correct, ctr, "ec: ec_add_test_crosscheck_double FAILED");
}

void ec_add_test_order_of_subgroup_is_order_of_generator(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X51441C4EE272FE55, 0X2DB9775DAEDDE550), bf_create_elem(0X12DD1A65F1D5B480, 0X6CB9034E20AD0EEB)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7AA01DC08C73455A, 0X51F2DF8B2F5FA18C), bf_create_elem(0X6730BC49B9A98F41, 0X57DEBE6DBCE321DE)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X80, 0), bf_create_elem(0X0000000200000000, 0X2000000000)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //(order-1) * GEN

	//Act
	ec_point_lproj sum = ec_add((ec_point_lproj) GEN, P);

	//Assert
	uint64_t correct = ec_equal_point_lproj((ec_point_lproj) INFTY, sum);
	assert_true(correct, ctr, "ec: ec_add_test_order_of_subgroup_is_order_of_generator FAILED");
}

void ec_double_test_example(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //78632917462800214 * GEN

	ef_intrl_elem EX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XE232234BCDA1F7A9, 0X72E20184752C18C3), bf_create_elem(0XA249666BD031EF41,0X78BFF5D2CAFC6FC2)));
	ef_intrl_elem EL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X1064B65732340CD9, 0X63C8BCD9FBA1773F), bf_create_elem(0X128743886A3EC579,0X6E8E5EC1D3091E7F)));
	ef_intrl_elem EZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X212FB30BAC886248, 0X797A835808F9BF57), bf_create_elem(0X8CA3223B872913CE,0X4B0C008D15FE6EB1)));
	ec_point_lproj expected = ec_create_point_lproj(EX, EL, EZ);

	//Act
	ec_point_lproj actual = ec_double(P);

	//Assert
	uint64_t correct = ec_equal_point_lproj(expected, actual) && ec_is_on_curve(P) && ec_is_on_curve(actual);
	assert_true(correct, ctr, "ec: ec_double_test_example FAILED");
}

void ec_double_test_is_on_curve_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();

		//Act
		ec_point_lproj P_doubled = ec_double(P);

		//Assert
		correct = ec_is_on_curve(P_doubled);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_double_test_is_on_curve_rnd FAILED");
}

void ec_double_test_double_of_sum_is_sum_of_doubled(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //78632917462800214 * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //99921481365893197563 * GEN

	//Act
	ec_point_lproj double_of_sum = ec_double(ec_add(P, Q));
	ec_point_lproj sum_of_doubled = ec_add(ec_double(P), ec_double(Q));

	//Assert
	uint64_t correct = ec_equal_point_lproj(double_of_sum, sum_of_doubled) && ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(double_of_sum);
	assert_true(correct, ctr, "ec: ec_double_test_double_of_sum_is_sum_of_doubled FAILED");
}

void ec_double_test_double_of_sum_is_sum_of_doubled_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj Q = ec_rand_point_lproj();

		//Act
		ec_point_lproj sum = ec_add(P, Q);
		ec_point_lproj double_of_sum = ec_double(sum);
		ec_point_lproj P_doubled = ec_double(P);
		ec_point_lproj Q_doubled = ec_double(Q);
		ec_point_lproj sum_of_doubled = ec_add(P_doubled, Q_doubled);

		//Assert
		correct = ec_equal_point_lproj(double_of_sum, sum_of_doubled);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			printf("Q: ");
			ec_print_hex(Q);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_double_test_double_of_sum_is_sum_of_doubled_rnd FAILED");
}

void ec_double_test_double_of_infty_is_infty(test_ctr *ctr) {
	//Arrange & Act
	ec_point_lproj doubled = ec_double((ec_point_lproj) INFTY);

	//Assert
	uint64_t correct = ec_equal_point_lproj((ec_point_lproj) INFTY, doubled);
	assert_true(correct, ctr, "ec: ec_double_test_double_of_infty_is_infty FAILED");
}

void ec_neg_test_example(test_ctr* ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB55AF6853BAA916A, 0X368B1A5434DF7331), bf_create_elem(0X4B3628231E3A83C2, 0XCA5DED000C90027)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //2243156791409331652485 * P

	ef_intrl_elem EX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4)));
	ef_intrl_elem EL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XB55AF6853BAA91FA, 0X368B1A5434DF7331), bf_create_elem(0X4B36A8231E3A83C2, 0XCA5DED000C90427)));
	ef_intrl_elem EZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400)));
	ec_point_lproj expected = ec_create_point_lproj(EX, EL, EZ);

	//Act
	ec_point_lproj actual = ec_neg(P);

	//Assert
	uint64_t correct = ec_equal_point_lproj(expected, actual) && ec_is_on_curve(P) && ec_is_on_curve(actual);
	assert_true(correct, ctr, "ec: ec_neg_test_example FAILED");
}

void ec_neg_test_is_on_curve_rnd(test_ctr* ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();

		//Act
		ec_point_lproj P_inv = ec_neg(P);

		//Assert
		correct = ec_is_on_curve(P_inv);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_neg_test_is_on_curve_rnd FAILED");
}

void ec_neg_test_neg_of_sum_same_as_sum_of_neg(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X51441C4EE272FE55, 0X2DB9775DAEDDE550), bf_create_elem(0X12DD1A65F1D5B480, 0X6CB9034E20AD0EEB)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7AA01DC08C73455A, 0X51F2DF8B2F5FA18C), bf_create_elem(0X6730BC49B9A98F41, 0X57DEBE6DBCE321DE)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X80, 0), bf_create_elem(0X0000000200000000, 0X2000000000)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //(order-1) * GEN

	ef_intrl_elem QX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266)));
	ef_intrl_elem QL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490)));
	ef_intrl_elem QZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001)));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //99921481365893197563 * GEN

	//Act
	ec_point_lproj neg_of_sum = ec_neg(ec_add(P, Q));
	ec_point_lproj sum_of_neg = ec_add(ec_neg(P), ec_neg(Q));

	//Assert
	uint64_t correct = ec_equal_point_lproj(neg_of_sum, sum_of_neg);
	assert_true(correct, ctr, "ec: ec_neg_test_neg_of_sum_same_as_sum_of_neg FAILED");
}

void ec_neg_test_neg_of_sum_same_as_sum_of_neg_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj Q = ec_rand_point_lproj();

		//Act
		ec_point_lproj neg_of_sum = ec_neg(ec_add(P, Q));
		ec_point_lproj sum_of_neg = ec_add(ec_neg(P), ec_neg(Q));

		//Assert
		correct = ec_equal_point_lproj(neg_of_sum, sum_of_neg);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			printf("Q: ");
			ec_print_hex(Q);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_neg_test_neg_of_sum_same_as_sum_of_neg_rnd FAILED");
}

void ec_neg_test_neg_of_double_same_as_double_of_neg(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //99921481365893197563 * GEN

	//Act
	ec_point_lproj neg_of_double = ec_neg(ec_double(P));
	ec_point_lproj double_of_neg = ec_double(ec_neg(P));

	//Assert
	uint64_t correct = ec_equal_point_lproj(neg_of_double, double_of_neg);
	assert_true(correct, ctr, "ec: ec_neg_test_neg_of_double_same_as_double_of_neg FAILED");
}

void ec_neg_test_neg_of_double_same_as_double_of_neg_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();

		//Act
		ec_point_lproj neg_of_double = ec_neg(ec_double(P));
		ec_point_lproj double_of_neg = ec_double(ec_neg(P));

		//Assert
		correct = ec_equal_point_lproj(neg_of_double, double_of_neg);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_neg_test_neg_of_double_same_as_double_of_neg_rnd FAILED");
}

void ec_neg_test_neg_of_neg_is_original(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //99921481365893197563 * GEN

	//Act
	ec_point_lproj P_neg = ec_neg(P);
	ec_point_lproj P_neg_neg = ec_neg(P_neg);

	//Assert
	uint64_t correct = ec_equal_point_lproj(P, P_neg_neg) && ec_is_on_curve(P) && ec_is_on_curve(P_neg);
	assert_true(correct, ctr, "ec: ec_neg_test_neg_of_neg_is_original FAILED");
}

void ec_neg_test_neg_of_neg_is_original_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();

		//Act
		ec_point_lproj P_neg = ec_neg(P);
		ec_point_lproj P_neg_neg = ec_neg(P_neg);

		//Assert
		correct = ec_equal_point_lproj(P, P_neg_neg);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_neg_test_neg_of_neg_is_original_rnd FAILED");
}

void ec_neg_of_infty_is_infty(test_ctr *ctr) {
	//Arrange & Act
	ec_point_lproj neg = ec_neg((ec_point_lproj) INFTY);

	//Assert
	uint64_t correct = ec_equal_point_lproj((ec_point_lproj) INFTY, neg);
	assert_true(correct, ctr, "ec: ec_double_test_neg_of_infty_is_infty FAILED");
}

void ec_neg_test_add_point_with_neg_is_infty(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000)));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //78632917462800214 * GEN

	//Act
	ec_point_lproj P_neg = ec_neg(P);
	ec_point_lproj sum = ec_add(P, P_neg);

	//Assert
	uint64_t correct = ec_equal_point_lproj((ec_point_lproj) INFTY, sum);
	assert_true(correct, ctr, "ec: ec_neg_test_add_point_with_neg_is_infty FAILED");
}

void ec_neg_test_add_point_with_neg_is_infty_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange
		ec_point_lproj P = ec_rand_point_lproj();

		//Act
		ec_point_lproj P_neg = ec_neg(P);
		ec_point_lproj sum = ec_add(P, P_neg);

		//Assert
		correct = ec_equal_point_lproj(sum, (ec_point_lproj) INFTY);
		if(!correct) {
			printf("P: ");
			ec_print_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_neg_test_add_point_with_inverse_is_infty_rnd FAILED");
}

void ec_add_mixed_test_crosscheck_rnd(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};
	ec_point_laffine P = ec_rand_point_laffine();
	ec_point_lproj P_proj = ec_create_point_lproj(P.x, P.l, one);
	ec_point_lproj Q = ec_rand_point_lproj();

	//Act
	ec_point_lproj expected = ec_add(P_proj, Q);
	ec_point_lproj result = ec_add_mixed(P, Q);

	//Assert
	uint64_t equal = ec_equal_point_lproj(expected, result);
	assert_true(equal, ctr, "ec: ec_add_mixed_test_is_on_curve_rnd FAILED");
}

void ec_double_then_add_test_is_on_curve_rnd(test_ctr *ctr) {
	//Arrange
	ec_point_laffine P = ec_rand_point_laffine();
	ec_point_lproj Q = ec_rand_point_lproj();

	//Act
	ec_point_lproj result = ec_double_then_add(P, Q);

	//Assert
	uint64_t on_curve = ec_is_on_curve(result);
	assert_true(on_curve, ctr, "ec: ec_double_then_add_test_is_on_curve_rnd FAILED");
}

void ec_double_then_addtwo_test_crosscheck_rnd(test_ctr *ctr) {
	//Arrange
	ef_intrl_elem one = (ef_intrl_elem) {{{1, 0}, {0, 0}}};
	ec_point_laffine P1 = ec_rand_point_laffine();
	ec_point_lproj P1proj = ec_create_point_lproj(P1.x, P1.l, one);
	ec_point_laffine P2 = ec_rand_point_laffine();
	ec_point_lproj P2proj = ec_create_point_lproj(P2.x, P2.l, one);
	ec_point_lproj Q = ec_rand_point_lproj();

	//Act
	ec_point_lproj expected = ec_add_mixed(P2, ec_double_then_add(P1, Q));
	
	ec_point_lproj result = ec_double_then_addtwo(P1, P2, Q);

	//Assert
	uint64_t equal = ec_equal_point_lproj(expected, result);
	assert_true(equal, ctr, "ec: ec_double_then_addtwo_test_is_on_curve_rnd FAILED");
}

void ec_endo_laffine_test_is_on_curve_rnd(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 3; i++) {
		//Arrange
		ec_point_laffine P = ec_rand_point_laffine();
		
		//Act
		ec_point_laffine phi_P = ec_endo_laffine(P);
		
		//Assert
		correct = ec_is_on_curve(ec_laffine_to_lproj(phi_P));
		if(!correct) {
			printf("P: \n");
			ec_print_hex_laffine(P);
			break;
		}
	}
	assert_true(correct, ctr, "ec: ec_endo_laffine_test_is_on_curve_rnd FAILED");
}

void ec_tests(test_ctr *ctr) {
	ec_create_point_lproj_test_example(ctr);

	ec_equal_point_lproj_test_equivalent_example(ctr);
	ec_equal_point_lproj_test_notequivalent_example(ctr);

	ec_is_on_curve_test_generator_on_curve(ctr);
	ec_is_on_curve_test_point_not_on_curve(ctr);
	ec_rand_point_lproj_test_on_curve(ctr);
	ec_rand_point_laffine_test_on_curve(ctr);

	ec_add_test_example(ctr);
	ec_add_test_is_on_curve_rnd(ctr);
	ec_add_test_associative(ctr);
	ec_add_test_associative_rnd(ctr);
	ec_add_test_commutative(ctr);
	ec_add_test_commutative_rnd(ctr);
	ec_add_test_point_plus_infty_is_point(ctr);
	ec_add_test_infty_plus_point_is_point(ctr);
	ec_add_test_crosscheck_double(ctr);
	ec_add_test_order_of_subgroup_is_order_of_generator(ctr);

	ec_double_test_example(ctr);
	ec_double_test_is_on_curve_rnd(ctr);
	ec_double_test_double_of_sum_is_sum_of_doubled(ctr);
	ec_double_test_double_of_sum_is_sum_of_doubled_rnd(ctr);
	ec_double_test_double_of_infty_is_infty(ctr);

	ec_neg_test_example(ctr);
	ec_neg_test_is_on_curve_rnd(ctr);
	ec_neg_test_neg_of_sum_same_as_sum_of_neg(ctr);
	ec_neg_test_neg_of_sum_same_as_sum_of_neg_rnd(ctr);
	ec_neg_test_neg_of_double_same_as_double_of_neg(ctr);
	ec_neg_test_neg_of_double_same_as_double_of_neg_rnd(ctr);
	ec_neg_test_neg_of_neg_is_original(ctr);
	ec_neg_test_neg_of_neg_is_original_rnd(ctr);
	ec_neg_of_infty_is_infty(ctr);
	ec_neg_test_add_point_with_neg_is_infty(ctr);
	ec_neg_test_add_point_with_neg_is_infty_rnd(ctr);

	ec_double_then_add_test_is_on_curve_rnd(ctr);
	 
	ec_add_mixed_test_crosscheck_rnd(ctr);
	
	ec_double_then_addtwo_test_crosscheck_rnd(ctr);
	
	ec_endo_laffine_test_is_on_curve_rnd(ctr);
}
