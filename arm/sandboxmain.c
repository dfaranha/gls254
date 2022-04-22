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
	
	int w = 5;
	int l = 33;
	
	uint64x2x2_t k = (uint64x2x2_t) SUBGROUP_ORDER;
	k.val[0][0]--;
	ec_split_scalar decomp = ec_scalar_decomp(k);
	decomp.k1[0]++;
	decomp.k2[0]++;
	
	signed char rec_k1[l];
	signed char rec_k2[l];
	reg_rec(decomp.k1, w, rec_k1, l-1);
	reg_rec(decomp.k2, w, rec_k2, l-1);
	
	printf("k1 sign: %lu\n", decomp.k1_sign);
	ec_print_rec(rec_k1, l);
	printf("k2 sign: %lu\n", decomp.k2_sign);
	ec_print_rec(rec_k2, l);
	
	dispose_components();
	return 0;
}
