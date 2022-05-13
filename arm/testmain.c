#include <stdio.h>

#include "tests/test_tool.h"

#include "common/setup.h"
#include "tests/basefield_tests.h"
#include "tests/benchmark_tool_tests.h"
#include "tests/ec_tests.h"
#include "tests/ec_scalarmull_tests.h"
#include "tests/extensionfield_tests.h"
#include "tests/extensionfield_interleaved_tests.h"
#include "tests/util_tests.h"

int main() {
	init_components();
	
	test_ctr *ctr = init_ctr();

	benchmark_tool_tests(ctr);
	util_tests(ctr);
	basefield_tests(ctr);
	extensionfield_tests(ctr);
	extensionfield_interleaved_tests(ctr);
	ec_tests(ctr);
	ec_scalarmull_tests(ctr);

	print_results(ctr);
	
	dispose_components();
	return 0;
}
