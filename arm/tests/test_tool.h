#include <arm_neon.h>

#ifndef TEST_TOOL_H
#define TEST_TOOL_H

typedef struct test_ctr_st {
	uint64_t num_run, num_passed;
} test_ctr;

test_ctr *init_ctr();

void assert_true(uint64_t cond, test_ctr *ctr, char * err_msg);

void assert_false(uint64_t cond, test_ctr *ctr, char * err_msg);

void print_results(test_ctr *ctr);

#endif
