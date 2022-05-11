#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/ec.h"
#include "common/ec_scalarmull.h"
#include "common/setup.h"
#include "common/utils.h"

int main() {
	init_components();

	ec_point_laffine P = ec_rand_point_laffine();
	uint64x2x2_t k = ec_rand_scalar();
	ec_point_laffine R;
	ec_scalarmull_single_endo_w4_table2D_bulk_ptr(&P, k, &R);
	ec_print_hex_laffine(R);
	
	dispose_components();
	return 0;
}
