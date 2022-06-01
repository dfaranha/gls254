#include <stdio.h>

#include "../common/extensionfield_interleaved.h"
#include "../common/utils.h"
#include "benchmark_extensionfield_interleaved.h"
#include "benchmark_tool.h"

void benchmark_ef_intrl_add() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c = ef_intrl_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_add_ptr() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c;
		ef_intrl_add_ptr(&a, &b, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_add_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_square() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c = ef_intrl_square(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_square\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_square_ptr() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c;
		ef_intrl_square_ptr(&a, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_square_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_mull() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c = ef_intrl_mull(a,b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_mull\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_mull_ptr() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		ef_intrl_elem b = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c;
		ef_intrl_mull_ptr(&a, &b, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_mull_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_red() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem_unred a = ef_intrl_rand_unred_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c = ef_intrl_red(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_red\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_red_ptr() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem_unred a = ef_intrl_rand_unred_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c;
		ef_intrl_red_ptr(&a, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_red_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c = ef_intrl_inv(a, 0);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_inv not in const time\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_inv_in_const() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_elem c = ef_intrl_inv(a, 1);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, c);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_inv in const time\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_sim_inv(uint64_t len) {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs];
	ef_intrl_elem sum = ef_intrl_rand_elem();
	ef_intrl_elem inputs[len];
	ef_intrl_elem outputs[len];

	for(int i = 0; i < num_runs; i++) {
		ef_intrl_elem a = ef_intrl_rand_elem();
		uint64_t start = read_pmccntr();
		ef_intrl_sim_inv(inputs, outputs, len, 0);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_intrl_add(sum, outputs[i % len]);
	}
	ef_intrl_print_hex_nl(sum);
	printf("BENCHMARK ef_intrl_sim_inv %lu elements\n", len);
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_intrl_all() {
	benchmark_ef_intrl_add();
	benchmark_ef_intrl_add_ptr();
	benchmark_ef_intrl_square();
	benchmark_ef_intrl_square_ptr();
	benchmark_ef_intrl_mull();
	benchmark_ef_intrl_mull_ptr();
	benchmark_ef_intrl_red();
	benchmark_ef_intrl_red_ptr();
	benchmark_ef_intrl_inv();
	benchmark_ef_intrl_inv_in_const();
	benchmark_ef_intrl_sim_inv(16);
	benchmark_ef_intrl_sim_inv(4);
}
