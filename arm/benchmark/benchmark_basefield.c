#include <stdio.h>

#include "../common/basefield.h"
#include "../common/utils.h"
#include "benchmark_basefield.h"
#include "benchmark_tool.h"

void benchmark_bf_add() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_pmull() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2x2_t c = bf_pmull(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, bf_red(c));
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_pmull\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_psquare() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2x2_t c = bf_psquare(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, bf_red_psquare(c));
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_psquare\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_red() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		poly64x2x2_t c = bf_pmull(a,b); //To get more avrg input
		uint64_t start = read_pmccntr();
		poly64x2_t d = bf_red(c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, d);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_red\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_red_lazy() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		poly64x2x2_t c = bf_pmull(a,b);
		uint64_t start = read_pmccntr();
		poly64x2_t d = bf_red_lazy(c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, d);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_red_lazy\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_red_psquare() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2x2_t c = bf_psquare(a); //To get more avrg input
		uint64_t start = read_pmccntr();
		poly64x2_t d = bf_red_psquare(c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, d);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_red_psquare\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_square() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t d = bf_square(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, d);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_square\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_mull() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t d = bf_mull(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, d);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_mull\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_addchain_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_addchain_inv(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_addchain_inv\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_addchain_lookup_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	precomp_inv_tables();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_addchain_lookup_inv(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_addchain_lookup_inv\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_loop_6() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_loop(a, 6);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_loop_6\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_lookup_6() {
	uint64_t num_runs = 200;
	uint64_t times[num_runs]; 
	precomp_inv_table(6);
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_lookup_6(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_lookup_6\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_loop_12() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_loop(a, 12);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_loop_12\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_lookup_12() {
	uint64_t num_runs = 200;
	uint64_t times[num_runs]; 
	precomp_inv_table(12);
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_lookup_12(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_lookup_12\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_loop_18() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_loop(a, 18);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_loop_18\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_lookup_18() {
	uint64_t num_runs = 200;
	uint64_t times[num_runs]; 
	precomp_inv_table(18);
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_lookup_18(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_lookup_18\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_loop_30() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_loop(a, 30);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_loop_30\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_lookup_30() {
	uint64_t num_runs = 200;
	uint64_t times[num_runs]; 
	precomp_inv_table(30);
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_lookup_30(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_lookup_30\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_loop_48() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_loop(a, 48);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_loop_48\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_multisquare_lookup_48() {
	uint64_t num_runs = 200;
	uint64_t times[num_runs]; 
	precomp_inv_table(48);
	poly64x2_t sum = bf_rand_elem();
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		poly64x2_t c = bf_multisquare_lookup_48(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = bf_add(sum, c);
	}
	bf_print_hex_nl(sum);
	printf("BENCHMARK bf_multisquare_lookup_48\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_all() {
	benchmark_bf_add();
	benchmark_bf_pmull();
	benchmark_bf_psquare();
	benchmark_bf_red();
	benchmark_bf_red_psquare();
	benchmark_bf_red_lazy();
	benchmark_bf_square();
	benchmark_bf_mull();
	benchmark_bf_addchain_inv();
	benchmark_bf_addchain_lookup_inv();
	/*benchmark_bf_multisquare_loop_6();
	benchmark_bf_multisquare_lookup_6();
	benchmark_bf_multisquare_loop_12();
	benchmark_bf_multisquare_lookup_12();
	benchmark_bf_multisquare_loop_18();
	benchmark_bf_multisquare_lookup_18();
	benchmark_bf_multisquare_loop_30();
	benchmark_bf_multisquare_lookup_30();
	benchmark_bf_multisquare_loop_48();
	benchmark_bf_multisquare_lookup_48();*/
}

