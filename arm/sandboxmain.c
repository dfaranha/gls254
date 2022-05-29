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


	//uint64x2x2_t k = (uint64x2x2_t) {{{1984, 0}, {0, 0}}};
	uint64x2x2_t k = (uint64x2x2_t) {{{4376904538686917875, 4547025796284947885}, {0, 2305843009213693952}}}; //r-2

	ec_split_scalar res = ec_scalar_decomp(k);
	printf("k1 sign = %lu [0] = %lu [1] = %lu\n", res.k1_sign, res.k1[0], res.k1[1]);
	printf("k2 sign = %lu [0] = %lu [1] = %lu\n", res.k2_sign, res.k2[0], res.k2[1]);

	int l = 64;
	int w = 3;
	uint64x2_t k1;
	k1[0] = 11269849258554726761U;
	k1[1] = 4611686018427387904U;
	signed char rec_k1[l];
	reg_rec(k1, w, rec_k1, l);

	dispose_components();
	return 0;
}
