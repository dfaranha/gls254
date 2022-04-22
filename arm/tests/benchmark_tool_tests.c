#include "../benchmark/benchmark_tool.h"
#include "benchmark_tool_tests.h"
#include <stdio.h>

void insert_sorted_test(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {2, 4, 7, 9, 0, 0};
	
	//Act
	insert_sorted(1, times, 4);
	
	//Assert
	uint64_t correctly_inserted = 
		(times[0] == 1) &&
		(times[1] == 2) &&
		(times[2] == 4) &&
		(times[3] == 7) &&
		(times[4] == 9) &&
		(times[5] == 0);
	assert_true(correctly_inserted, ctr, "benchmark_tool: insert_sorted_test FAILED");
}

void insert_sorted_test_notzeroinit(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {2, 4, 7, 9, 10, 1};
	
	//Act
	insert_sorted(1, times, 4);
	
	//Assert
	uint64_t correctly_inserted = 
		(times[0] == 1) &&
		(times[1] == 2) &&
		(times[2] == 4) &&
		(times[3] == 7) &&
		(times[4] == 9) &&
		(times[5] == 1);
	assert_true(correctly_inserted, ctr, "benchmark_tool: insert_sorted_test_notzeroinit FAILED");
}

void insert_sorted_test_len1(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {3};
	
	//Act
	insert_sorted(4, times, 0);
	
	//Assert
	uint64_t correctly_inserted = times[0] == 4;
	assert_true(correctly_inserted, ctr, "benchmark_tool: insert_sorted_test_len1 FAILED");
}

void benchmark_tool_tests(test_ctr *ctr) {
	insert_sorted_test(ctr);
	insert_sorted_test_notzeroinit(ctr);
	insert_sorted_test_len1(ctr);
}
