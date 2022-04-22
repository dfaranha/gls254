#include <arm_neon.h>
#include <time.h>

#ifndef BENCHMARK_TOOL_H
#define BENCHMARK_TOOL_H

uint64_t read_pmccntr();

int64_t cpu_nseconds(void);

void insert_sorted(uint64_t time, uint64_t *sorted_times, uint64_t curr_len);

#endif
