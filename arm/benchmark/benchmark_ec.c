#include <stdio.h>

#include "../common/ec.h"
#include "../common/utils.h"
#include "benchmark_ec.h"
#include "benchmark_tool.h"


void benchmark_ec_add() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj a = ec_rand_point_lproj();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_unchecked() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj a = ec_rand_point_lproj();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_add_unchecked(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_unchecked\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_mixed() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_add_mixed(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_mixed\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_mixed_unchecked() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_add_mixed_unchecked(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_mixed_unchecked\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj a = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_double(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_mixed() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_double_mixed(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_mixed\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_alt() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj a = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_double_alt(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_alt\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_then_add() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_double_then_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_then_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_then_addtwo() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_laffine b = ec_rand_point_laffine();
		ec_point_lproj c = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj d = ec_double_then_addtwo(a, b, c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, d);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_then_addtwo\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_then_add_nonatomic() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_double_then_add_nonatomic(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_then_add_nonatomic\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_then_addtwo_nonatomic() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_laffine b = ec_rand_point_laffine();
		ec_point_lproj c = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj d = ec_double_then_addtwo_nonatomic(a, b, c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, d);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_then_addtwo_nonatomic\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_all() {
	benchmark_ec_add();
	benchmark_ec_add_unchecked();
	benchmark_ec_add_mixed();
	benchmark_ec_add_mixed_unchecked();
	benchmark_ec_double();
	benchmark_ec_double_mixed();
	benchmark_ec_double_alt();
	benchmark_ec_double_then_add();
	benchmark_ec_double_then_addtwo();
	benchmark_ec_double_then_add_nonatomic();
	benchmark_ec_double_then_addtwo_nonatomic();
}
