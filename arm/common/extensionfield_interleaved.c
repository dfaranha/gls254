#include "extensionfield_interleaved.h"

#include "basefield.h"
#include "extensionfield.h"
#include "utils.h"
#include <stdio.h>

/* Field elements are now shuffled like this:
 * a = {{a[0][0], a[1][0]}, {a[0][1], a[1][1]}}
 */

void ef_intrl_print_expr(ef_intrl_elem a) {
	ef_print_expr(ef_intrl_disentangle(a));
}

void ef_intrl_print_expr_nl(ef_intrl_elem a) {
	ef_intrl_print_expr(a);
	printf("\n");
}

void ef_intrl_print_hex(ef_intrl_elem a) {
	ef_print_hex(ef_intrl_disentangle(a));
}

void ef_intrl_print_hex_nl(ef_intrl_elem a) {
	ef_intrl_print_hex(a);
	printf("\n");
}

void ef_intrl_print_unred_expr(ef_intrl_elem_unred a) {
	poly64x2x2_t a0 = ef_intrl_disentangle_unred_lower(a);
	poly64x2x2_t a1 = ef_intrl_disentangle_unred_higher(a);
	printf("(");
	bf_print_unred_expr(a1);
	printf(")u + (");
	bf_print_unred_expr(a0);
	printf(")");
}

void ef_intrl_print_unred_expr_nl(ef_intrl_elem_unred a) {
	ef_intrl_print_unred_expr(a);
	printf("\n");
}

ef_elem ef_intrl_disentangle(ef_intrl_elem a) {
	ef_elem res;
	poly64x2_t t = vextq_p64(a.val[0], a.val[0], 1); //t = {a[1][0], a[0][0]}
	res.val[0] = vextq_p64(t, a.val[1], 1); //{a[0][0], a[0][1]}
	res.val[1] = a.val[1]; //{a[0][1], a[1][1]}
	res.val[1][0] = t[0]; //{a[1][0], a[1][1]}
	return res;
}

//{{a[0][0], a[0][1]}, {a[1][0], a[1][1]}}
ef_intrl_elem ef_intrl_interleave(ef_elem a) {
	ef_intrl_elem res;
	poly64x2_t t = vextq_p64(a.val[0], a.val[0], 1); //t = {a[0][1], a[0][0]}
	res.val[0] = vextq_p64(t, a.val[1], 1); //{a[0][0], a[1][0]}
	res.val[1] = a.val[1]; //{a[1][0], a[1][1]}
	res.val[1][0] = t[0]; //{a[0][1], a[1][1]}
	return res;
}

poly64x2x2_t ef_intrl_disentangle_unred_lower(ef_intrl_elem_unred c) {
	return (poly64x2x2_t) {{{c.val[0][0], c.val[1][0]}, {c.val[2][0], c.val[3][0]}}};
}

poly64x2x2_t ef_intrl_disentangle_unred_higher(ef_intrl_elem_unred c) {
	return (poly64x2x2_t) {{{c.val[0][1], c.val[1][1]}, {c.val[2][1], c.val[3][1]}}};
}

ef_intrl_elem ef_intrl_rand_elem() {
	return ef_intrl_interleave(ef_rand_elem());
}

ef_intrl_elem_unred ef_intrl_rand_unred_elem() {
	ef_intrl_elem_unred c;

	c.val[0] = bf_rand_elem();
	c.val[1] = bf_rand_elem();
	c.val[2] = bf_rand_elem();
	c.val[3] = bf_rand_elem();

	return c;
}

ef_intrl_elem ef_intrl_red_from_lazy(ef_intrl_elem c) {
	ef_intrl_elem res;
	poly64x2_t mask = vdupq_n_p64(pow2to63);
	poly64x2_t bit127 = (poly64x2_t) vandq_u64((uint64x2_t) c.val[1], (uint64x2_t) mask);
	res.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[1], (uint64x2_t) bit127); //Set bit127 to 0
	res.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) c.val[0], (uint64x2_t) bit127); //Add z^63 if bit127 set
	bit127 = (poly64x2_t) vshrq_n_u64((uint64x2_t) bit127, 63);
	res.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) res.val[0], (uint64x2_t) bit127); //Add 1 if bit127 set
	return res;
}

uint64_t ef_intrl_equal(ef_intrl_elem a, ef_intrl_elem b) {
	ef_intrl_elem a_red = ef_intrl_red_from_lazy(a);
	ef_intrl_elem b_red = ef_intrl_red_from_lazy(b);
	return equal_poly64x2x2(a_red, b_red);
}

ef_intrl_elem ef_intrl_inv(ef_intrl_elem a) {
	ef_intrl_elem_unred c;
	poly64x2x2_t r0, r1;
	poly64x2_t t, t0, t1, t2, t3, t4, t5;
	poly64x2_t z = vdupq_n_p64(0);

	t2 = vextq_p64(a.val[0], a.val[0], 1); //{a1[0], a0[0]}
	t3 = vextq_p64(a.val[1], a.val[1], 1); //{a1[1], a0[1]}

	//a0 + a1
	t4 = (poly64x2_t) veorq_u64((uint64x2_t) t2, (uint64x2_t) a.val[0]); //{a0[0] ^ a1[0], a0[0] ^ a1[0]}
	t5 = (poly64x2_t) veorq_u64((uint64x2_t) t3, (uint64x2_t) a.val[1]); //{a0[1] ^ a1[1], a0[1] ^ a1[1]}

	//a0^2 + a1^2
	r0.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t4[0], t4[0]));
	r0.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t5[0], t5[0]));

	//a0 * a1
	r1.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], t2[0]));
	r1.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], t3[0]));
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[0][0], t3[0]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a.val[1][0], t2[0]));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	r1.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) r1.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	r1.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) r1.val[1], (uint64x2_t) t1);

	//t = a0 * a1 + a0^2 + a1^2
	r0.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) r0.val[0], (uint64x2_t) r1.val[0]);
	r0.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) r0.val[1], (uint64x2_t) r1.val[1]);

	//t^-1
	t = bf_inv(bf_red_lazy(r0));

	//c1 = a1*t^-1
	r1.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t[0], t2[0])); //t[0] * a1[0]
	r1.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(t, a.val[1])); //t[1] * a1[1]
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t[0], t3[0])); //t[0] * a1[1]
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(t, a.val[0])); //t[1] * a1[0]
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	r1.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) r1.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	r1.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) r1.val[1], (uint64x2_t) t1);

	//c0 = (a0 + a1)*t^-1
	r0.val[0] = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t[0], t4[0]));
	r0.val[1] = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(t, t5));
	t1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t[0], t5[0]));
	t0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_high_p64(t, t4));
	t0 = (poly64x2_t) veorq_u64((uint64x2_t) t0, (uint64x2_t) t1);
	t1 = vextq_p64(z, t0, 1);
	r0.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) r0.val[0], (uint64x2_t) t1);
	t1 = vextq_p64(t0, z, 1);
	r0.val[1] = (poly64x2_t) veorq_u64((uint64x2_t) r0.val[1], (uint64x2_t) t1);

	//Combining
	c.val[0] = (poly64x2_t) vzip1q_u64((uint64x2_t) r0.val[0], (uint64x2_t) r1.val[0]);
	c.val[1] = (poly64x2_t) vzip2q_u64((uint64x2_t) r0.val[0], (uint64x2_t) r1.val[0]);
	c.val[2] = (poly64x2_t) vzip1q_u64((uint64x2_t) r0.val[1], (uint64x2_t) r1.val[1]);
	c.val[3] = (poly64x2_t) vzip2q_u64((uint64x2_t) r0.val[1], (uint64x2_t) r1.val[1]);

	return ef_intrl_red(c);
}

void ef_intrl_sim_inv(ef_intrl_elem inputs[], ef_intrl_elem outputs[], uint64_t len) {
	ef_intrl_elem c[len];
	c[0] = inputs[0];
	for(int i = 1; i < len; i++) {
		c[i] = ef_intrl_mull(c[i-1], inputs[i]);
	}
	ef_intrl_elem u = ef_intrl_inv(c[len-1]); //(a0 * a1 * ... * ak)^-1
	for(int i = len-1; i >= 1; i--) {
		outputs[i] = ef_intrl_mull(u, c[i-1]); //(a0 * a1 * ... ai-1 * ai)^-1 * (a0 * a1 * ... * ai-1) = ai^-1
		u = ef_intrl_mull(u, inputs[i]);
	}
	outputs[0] = u;
}
