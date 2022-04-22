#include <stdio.h>
#include <stdlib.h>

#include "test_tool.h"

test_ctr *init_ctr() {
	test_ctr *ctr = (test_ctr *) malloc(sizeof(test_ctr));
	ctr->num_run = 0;
	ctr->num_passed = 0;
	return ctr;
}

void assert_true(uint64_t cond, test_ctr *ctr, char * err_msg) {
	ctr->num_run++;
	if (cond) {
		ctr->num_passed++;
	} else {
		printf("%s\n", err_msg);
	}
}

void assert_false(uint64_t cond, test_ctr *ctr, char * err_msg) {
	ctr->num_run++;
	if (cond) {
		printf("%s\n", err_msg);
	} else {
		ctr->num_passed++;
	}
}

void print_results(test_ctr *ctr) {
	printf("Passed %lu of %lu tests.\n", ctr->num_passed, ctr->num_run);
	free(ctr);
}




