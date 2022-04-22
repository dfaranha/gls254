#include <unistd.h>

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
	// benchmark_bf_all();
	// benchmark_ef_all();
	// benchmark_ef_intrl_all();
	//benchmark_ec_all();
	benchmark_ec_scalarmull_all();

	dispose_components();
	return 0;
}
