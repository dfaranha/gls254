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

	ef_intrl_elem PX = ef_intrl_interleave(ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U)));
	ef_intrl_elem PL = ef_intrl_interleave(ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U)));
	ef_intrl_elem PZ = ef_intrl_interleave(ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)));
	ec_point_laffine P = ec_lproj_to_laffine(ec_create_point_lproj(PX, PL, PZ)); //P = 12345 * GEN

	ec_point_laffine table[16];
	ec_precompute_w4_table2D_ptr(&P, table);

	ef_intrl_elem zero = (ef_intrl_elem) {{{0, 0}, {0, 0}}};
	ef_intrl_elem a = (ef_intrl_elem) {{{pow2to63 + 1, pow2to63 + 1}, {pow2to63 + 1, pow2to63 + 1}}};
	ec_point_laffine P1 = ec_create_point_laffine(zero, zero);
	ec_point_laffine P2 = ec_create_point_laffine(zero, zero);
	uint64_t index1 = 15;
	uint64_t index2 = 1;
	bf_print_hex_nl(table[index1].x.val[0]);
	bf_print_hex_nl(table[index1].x.val[1]);
	bf_print_hex_nl(table[index1].l.val[0]);
	bf_print_hex_nl(table[index1].l.val[1]);
	printf("\n");
	bf_print_hex_nl(table[index2].x.val[0]);
	bf_print_hex_nl(table[index2].x.val[1]);
	bf_print_hex_nl(table[index2].l.val[0]);
	bf_print_hex_nl(table[index2].l.val[1]);
	printf("\n");
	printf("\n");

	lin_pass_w4_table2D_bulk_neon(&P1, &P2, table, index1, index2);
	bf_print_hex_nl(P1.x.val[0]);
	bf_print_hex_nl(P1.x.val[1]);
	bf_print_hex_nl(P1.l.val[0]);
	bf_print_hex_nl(P1.l.val[1]);

	printf("\n");
	bf_print_hex_nl(P2.x.val[0]);
	bf_print_hex_nl(P2.x.val[1]);
	bf_print_hex_nl(P2.l.val[0]);
	bf_print_hex_nl(P2.l.val[1]);
	
	//ec_print_hex_laffine(P1);
	//ec_print_hex_laffine(P2);
	
	dispose_components();
	return 0;
}
