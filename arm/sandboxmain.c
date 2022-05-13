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

	ef_intrl_elem EX = ef_intrl_interleave(ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25)));
	ef_intrl_elem EL = ef_intrl_interleave(ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB)));
	ef_intrl_elem EZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)));
	ec_point_laffine P = ec_lproj_to_laffine(ec_create_point_lproj(EX, EL, EZ)); // expected = 1984 * GEN
	uint64x2x2_t k = (uint64x2x2_t) {{{1984, 0}, {0, 0}}};
	ec_point_laffine R;
	//printf("New\n");
	ec_scalarmull_single_endo_w4_table2D_ptr(&P, k, &R);
	//printf("\nOld\n");
	//ec_scalarmull_single_endo_w5_randaccess_ptr(&P, k, &R);
	
	dispose_components();
	return 0;
}
