#include <arm_neon.h>

#include "basefield.h"
#include "extensionfield.h"
#include "ec.h"
#include "ec_scalarmull.h"

#ifndef UTILS_H
#define UTILS_H

//c += a
#define ADDACC_128(a0, a1, c0, c1)\
	asm volatile("ADDS %0, %0, %2;"\
		 "ADC %1, %1, %3;"\
		 : "+r" (c0), "+r" (c1)\
		 : "r" (a0), "r" (a1)\
		 );

#define ADDACC_192(a0, a1, a2, c0, c1, c2)\
	asm volatile("ADDS %0, %0, %3;"\
		 "ADCS %1, %1, %4;"\
		 "ADC %2, %2, %5;"\
		 : "+r" (c0), "+r" (c1), "+r" (c2)\
		 : "r" (a0), "r" (a1), "r" (a2)\
		 );

#define ADDACC_256(a0, a1, a2, a3, c0, c1, c2, c3)\
	asm volatile("ADDS %0, %0, %4;"\
		 "ADCS %1, %1, %5;"\
		 "ADCS %2, %2, %6;"\
		 "ADC %3, %3, %7;"\
		 : "+r" (c0), "+r" (c1), "+r" (c2), "+r" (c3)\
		 : "r" (a0), "r" (a1), "r" (a2), "r" (a3)\
		 );

//c -= a
#define SUBACC_128(a0, a1, c0, c1)\
	asm volatile("SUBS %0, %0, %2;" \
		 "SBC %1, %1, %3;"\
		: "+r" (c0), "+r" (c1)\
		: "r" (a0), "r" (a1)\
		);

#define SUBACC_192(a0, a1, a2, c0, c1, c2)\
	asm volatile("SUBS %0, %0, %3;" \
		 "SBCS %1, %1, %4;"\
		 "SBC %2, %2, %5;"\
		: "+r" (c0), "+r" (c1), "+r" (c2)\
		: "r" (a0), "r" (a1), "r" (a2)\
		);

#define SUBACC_256(a0, a1, a2, a3, c0, c1, c2, c3)\
	asm volatile("SUBS %0, %0, %4;" \
		 "SBCS %1, %1, %5;"\
		 "SBCS %2, %2, %6;"\
		 "SBC %3, %3, %7;"\
		: "+r" (c0), "+r" (c1), "+r" (c2), "+r" (c3)\
		: "r" (a0), "r" (a1), "r" (a2), "r" (a3)\
		);

#define CMOV(tmpx1, cmpval, eqcondx1, old, new, old_ptrx1, new_ptrx1, type) \
	tmpx1[0] = cmpval & eqcondx1[0]; \
	tmpx1 = vceq_u64(tmpx1, eqcondx1); \
	old_ptrx1[0] = (uint64_t) &old; \
	new_ptrx1[0] = (uint64_t) &new; \
	tmpx1 = vbsl_u64(tmpx1, new_ptrx1, old_ptrx1); \
	old = *((type*) tmpx1[0]); \

#define CSEL(cmpval, eqcond, old, new, new_ptr, type) \
	asm volatile ("CMP %1, %2;" \
	"CSEL %0, %3, %4, EQ;" \
	: "+r" ( new_ptr ) \
	: "r" (cmpval), "r" (eqcond), "r" (&new), "r" (&old) \
	); \
	old = *((type*)new_ptr); \

void utils_init();

double average(uint64_t nums[], uint64_t len);

uint64_t compare_doubles(double a, double b, double errmargin);

uint64_t equal_poly64x2(poly64x2_t a, poly64x2_t b);

uint64_t equal_poly64x2x2(poly64x2x2_t a, poly64x2x2_t b);

poly64x2x2_t concat_bf_poly(poly64x2_t p0, poly64x2_t p1);

double median(uint64_t sorted_nums[], uint64_t len);

uint64_t rand_uint64();

ec_point_laffine lproj_to_laffine(ec_point_lproj P);

uint64x2_t mult_u64(uint64_t a, uint64_t b);

#endif
