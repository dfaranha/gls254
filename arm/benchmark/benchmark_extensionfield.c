#include <stdio.h>

#include "../common/extensionfield.h"
#include "../common/utils.h"
#include "benchmark_extensionfield.h"
#include "benchmark_tool.h"


void benchmark_ef_add() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	ef_elem sum = ef_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		ef_elem a = ef_rand_elem();
		ef_elem b = ef_rand_elem();
		uint64_t start = read_pmccntr();
		ef_elem c = ef_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_add(sum, c);
	}
	ef_print_hex_nl(sum);
	printf("BENCHMARK ef_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_mull() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	ef_elem sum = ef_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		ef_elem a = ef_rand_elem();
		ef_elem b = ef_rand_elem();
		uint64_t start = read_pmccntr();
		ef_elem c = ef_mull(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_add(sum, c);
	}
	ef_print_hex_nl(sum);
	printf("BENCHMARK ef_mull\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_mull_B() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	ef_elem sum = ef_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		ef_elem a = ef_rand_elem();
		uint64_t start = read_pmccntr();
		ef_elem c = ef_mull_B(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_add(sum, c);
	}
	ef_print_hex_nl(sum);
	printf("BENCHMARK ef_mull_B\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_square() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	ef_elem sum = ef_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		ef_elem a = ef_rand_elem();
		uint64_t start = read_pmccntr();
		ef_elem c = ef_square(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_add(sum, c);
	}
	ef_print_hex_nl(sum);
	printf("BENCHMARK ef_square\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	ef_elem sum = ef_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		ef_elem a = ef_rand_elem();
		uint64_t start = read_pmccntr();
		ef_elem c = ef_inv(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ef_add(sum, c);
	}
	ef_print_hex_nl(sum);
	printf("BENCHMARK ef_inv\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ef_all() {
	benchmark_ef_add();
	benchmark_ef_mull();
	benchmark_ef_mull_B();
	benchmark_ef_square();
	benchmark_ef_inv();
}
