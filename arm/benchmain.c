#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include "armpmu_lib.h"

#include "benchmark/benchmark_tool.h"
#include "benchmark/benchmark_basefield.h"
#include "benchmark/benchmark_extensionfield.h"
#include "benchmark/benchmark_extensionfield_interleaved.h"
#include "benchmark/benchmark_ec.h"
#include "benchmark/benchmark_ec_scalarmull.h"
#include "common/setup.h"

int main() {
	init_components();

	nice(-30);
	enable_pmu(0x008);
	//benchmark_bf_all();
	//benchmark_ef_all();
	//benchmark_ef_intrl_all();
	//benchmark_ec_all();
	benchmark_ec_scalarmull_all();
	disable_pmu(0x008);

	dispose_components();
	return 0;
}
