#include <stdio.h>
#include "benchmark_tool.h"
#include "../common/utils.h"

inline uint64_t read_pmccntr() {
		uint64_t val;
		asm volatile("mrs %0, pmccntr_el0" : "=r"(val));
		return val;
}

// https://github.com/microsoft/FourQlib/blob/ff61f680505c98c98e33387962223ce0b5e620bc/FourQ_ARM_side_channel/tests/test_extras.c#L19
int64_t cpu_nseconds(void) {
	struct timespec time;
	clock_gettime(CLOCK_REALTIME, &time);
	return (int64_t)(time.tv_sec*1e9 + time.tv_nsec);
}

void insert_sorted(uint64_t time, uint64_t *sorted_times, uint64_t curr_len) {
	for (int i = 0; i < curr_len; i++) {
		if (sorted_times[i] > time) {
			uint64_t tmp = time;
			time = sorted_times[i];
			sorted_times[i] = tmp;
		}
	}
	sorted_times[curr_len] = time;
}
