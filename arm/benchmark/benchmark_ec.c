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

void benchmark_ec_add_unchecked_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj Q = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj R;
		ec_add_unchecked_ptr(&P, &Q, &R);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, R);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_unchecked_ptr\n");
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

void benchmark_ec_add_mixed_unchecked_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c;
		ec_add_mixed_unchecked_ptr(&a, &b, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_mixed_unchecked_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_laffine_unchecked() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_laffine b = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj c = ec_add_laffine_unchecked(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_laffine_unchecked\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_laffine_unchecked_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_laffine b = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj c;
		ec_add_laffine_unchecked_ptr(&a, &b, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_laffine_unchecked_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_sub_laffine_unchecked_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		ec_point_laffine Q = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj Radd, Rsub;
		ec_add_sub_laffine_unchecked_ptr(&P, &Q, &Radd, &Rsub);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, Radd);
		sum = ec_add(sum, Rsub);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_sub_laffine_unchecked_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_sub_mixed_unchecked_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		ec_point_lproj Q = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj Radd, Rsub;
		ec_add_sub_mixed_unchecked_ptr(&P, &Q, &Radd, &Rsub);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, Radd);
		sum = ec_add(sum, Rsub);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_sub_mixed_unchecked_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_add_endo_laffine_unchecked_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj R;
		ec_add_endo_laffine_unchecked_ptr(&P, &R);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, R);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_add_endo_laffine_unchecked_ptr\n");
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

void benchmark_ec_double_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj a = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c; 
		ec_double_ptr(&a, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_ptr\n");
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

void benchmark_ec_double_mixed_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj c;
		ec_double_mixed_ptr(&a, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_mixed_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double_alt_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj P = ec_rand_point_lproj();
		ec_point_lproj R;
		uint64_t start = read_pmccntr();
		ec_double_alt_ptr(&P, &R);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, R);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_alt_ptr\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_triple_mixed_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		ec_point_lproj R;
		uint64_t start = read_pmccntr();
		ec_triple_mixed_ptr(&P, &R);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, R);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_triple_mixed_ptr\n");
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

void benchmark_ec_double_then_add_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_lproj b = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj c;
		ec_double_then_add_ptr(&a, &b, &c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, c);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_then_add_ptr\n");
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

void benchmark_ec_double_then_addtwo_ptr() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine a = ec_rand_point_laffine();
		ec_point_laffine b = ec_rand_point_laffine();
		ec_point_lproj c = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj d;
		ec_double_then_addtwo_ptr(&a, &b, &c, &d);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, d);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_double_then_addtwo_ptr\n");
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

void benchmark_ec_endo_laffine() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_laffine R = ec_endo_laffine(P);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_endo_laffine\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_endo_lproj() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = (ec_point_lproj) INFTY;

	for(int i = 0; i < num_runs; i++) {
		ec_point_lproj P = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj R = ec_endo_lproj(P);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, R);
	}
	ec_print_hex(sum);
	printf("BENCHMARK ec_endo_lproj\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_all() {
	//benchmark_ec_add();
	benchmark_ec_add_unchecked_ptr();
	//benchmark_ec_add_mixed();
	benchmark_ec_add_mixed_unchecked_ptr();
	benchmark_ec_add_laffine_unchecked_ptr();
	benchmark_ec_add_sub_laffine_unchecked_ptr();
	benchmark_ec_add_sub_mixed_unchecked_ptr();
	benchmark_ec_add_endo_laffine_unchecked_ptr();
	benchmark_ec_double_ptr();
	benchmark_ec_double_mixed_ptr();
	benchmark_ec_double_alt_ptr();
	benchmark_ec_triple_mixed_ptr();
	benchmark_ec_double_then_add_ptr();
	benchmark_ec_double_then_addtwo_ptr();
	//benchmark_ec_double_then_add_nonatomic();
	//benchmark_ec_double_then_addtwo_nonatomic();
	benchmark_ec_endo_laffine();
	benchmark_ec_endo_lproj();
}
