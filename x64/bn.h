typedef unsigned int uint128_t __attribute__((mode(TI)));

/* Schoolbook multiplication. */
#define RLC_MUL_DIG(H, L, A, B)												\
	H = ((uint128_t)(A) * (uint128_t)(B)) >> 64;							\
	L = (A) * (B);															\

#define RLC_COMBA_ADD(T, R2, R1, R0, A)										\
	(T) = (R1);																\
	(R0) += (A);															\
	(R1) += (R0) < (A);														\
	(R2) += (R1) < (T);														\

#define RLC_COMBA_STEP_MUL(R2, R1, R0, A, B)								\
	uint64_t _r, _r0, _r1;													\
	RLC_MUL_DIG(_r1, _r0, A, B);											\
	RLC_COMBA_ADD(_r, R2, R1, R0, _r0);										\
	(R1) += _r1;															\
	(R2) += (R1) < _r1;														\

uint64_t bn_add1_low(uint64_t *c, const uint64_t *a, uint64_t digit, int size) {
    int i;
	register uint64_t carry, r0;

	carry = digit;
	for (i = 0; i < size && carry; i++, a++, c++) {
		r0 = (*a) + carry;
		carry = (r0 < carry);
		(*c) = r0;
	}
	for (; i < size; i++, a++, c++) {
		(*c) = (*a);
	}
	return carry;
}

uint64_t bn_addn_low(uint64_t *c, const uint64_t *a, const uint64_t *b, int size) {
	register uint64_t carry, c0, c1, r0, r1;

	carry = 0;
	for (int i = 0; i < size; i++, a++, b++, c++) {
		r0 = (*a) + (*b);
		c0 = (r0 < (*a));
		r1 = r0 + carry;
		c1 = (r1 < r0);
		carry = c0 | c1;
		(*c) = r1;
	}
	return carry;
}

uint64_t bn_subn_low(uint64_t *c, const uint64_t *a, const uint64_t *b, int size) {
	register uint64_t carry, r0, diff;

	/* Zero the carry. */
	carry = 0;
	for (int i = 0; i < size; i++, a++, b++, c++) {
		diff = (*a) - (*b);
		r0 = diff - carry;
		carry = ((*a) < (*b)) || (carry && !diff);
		(*c) = r0;
	}
	return carry;
}

void bn_muln_low(uint64_t *c, const uint64_t *a, const uint64_t *b, int size) {
	int i, j;
	const uint64_t *tmpa, *tmpb;
	uint64_t r0, r1, r2;

	r0 = r1 = r2 = 0;
	for (i = 0; i < size; i++, c++) {
		tmpa = a;
		tmpb = b + i;
		for (j = 0; j <= i; j++, tmpa++, tmpb--) {
			RLC_COMBA_STEP_MUL(r2, r1, r0, *tmpa, *tmpb);
		}
		*c = r0;
		r0 = r1;
		r1 = r2;
		r2 = 0;
	}
	for (i = 0; i < size; i++, c++) {
		tmpa = a + i + 1;
		tmpb = b + (size - 1);
		for (j = 0; j < size - (i + 1); j++, tmpa++, tmpb--) {
			RLC_COMBA_STEP_MUL(r2, r1, r0, *tmpa, *tmpb);
		}
		*c = r0;
		r0 = r1;
		r1 = r2;
		r2 = 0;
	}
}

#include <inttypes.h>

void bn_copy_cond(uint64_t *c, const uint64_t *a, int digits, uint64_t cond) {
	uint64_t mask, t;

	mask = -cond;
	for (int i = 0; i < digits; i++) {
		t = (a[i] ^ c[i]) & mask;
		c[i] ^= t;
	}
}

void bn_print(uint64_t *a, int digits) {
	int i;

	/* Suppress possible unused parameter warning. */
	(void)a;
	for (i = digits - 1; i >= 0; i--) {
        printf("%.16" PRIX64, a[i]);
	}
	printf("\n");

	return;
}
