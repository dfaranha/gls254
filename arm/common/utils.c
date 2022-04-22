#include <stdlib.h>
#include <time.h>

#include "utils.h"

void utils_init() {
	// Intializes random number generator
	time_t t;
	srand((unsigned) time(&t));
}

double average(uint64_t nums[], uint64_t len) {
	double sum = 0.0;
	for (int i = 0; i < len; i++) {
		sum += nums[i];
	}
	return sum / len;
}

/* Returns:
 *  0 if diff within [-errmargin, errmargin],
 *  1 if a is greater
 * -1 if b is greater
 */
uint64_t compare_doubles(double a, double b, double errmargin) {
	double diff = a - b;
	if (diff > errmargin) {
		return 1;
	}
	if (diff < -errmargin) {
		return -1;
	}
	return 0;
}

uint64_t equal_poly64x2(poly64x2_t a, poly64x2_t b) {
	return (a[0] == b[0]) && (a[1] == b[1]);
}

uint64_t equal_poly64x2x2(poly64x2x2_t a, poly64x2x2_t b) {
	return (a.val[0][0] == b.val[0][0]) &&
		   (a.val[0][1] == b.val[0][1]) &&
		   (a.val[1][0] == b.val[1][0]) &&
		   (a.val[1][1] == b.val[1][1]);
}

poly64x2x2_t concat_bf_poly(poly64x2_t p0, poly64x2_t p1) {
	poly64x2x2_t p;
	p.val[0] = p0;
	p.val[1] = p1;
	return p;
}

double median(uint64_t sorted_nums[], uint64_t len) {
	if (len % 2 == 0) {
		return (sorted_nums[len / 2 - 1] + sorted_nums[len / 2]) / 2.0;
	}
	return sorted_nums[len / 2];
}

uint64_t rand_uint64() {
	uint64_t r = 0;
	for (int i=0; i<64; i++) {
		r = r*2 + rand()%2;
	}
	return r;
}

//We can't multiply 64 bit integers directly, so need this :/
//Read a paper that seemed to indicate Karatsuba is not rly worth it for small bitsizes.
//Inspiration from here: https://stackoverflow.com/questions/25095741/how-can-i-multiply-64-bit-operands-and-get-128-bit-result-portably
uint64x2_t mult_u64(uint64_t a, uint64_t b) {
	uint64x2_t result;
	uint64_t al = a & 0xFFFFFFFF;
	uint64_t ah = a >> 32;
	uint64_t bl = b & 0XFFFFFFFF;
	uint64_t bh = b >> 32;
	uint64_t albl = al * bl;
	uint64_t albh = al * bh;
	uint64_t ahbl = ah * bl;
	uint64_t ahbh = ah * bh;
	
	//Can see that this way we don't get carry that must be handled :D
	uint64_t cross = (albl >> 32) + (ahbl & 0xFFFFFFFF) + albh;
	result[0] = (cross << 32) | (albl & 0xFFFFFFFF);
    result[1] = (ahbl >> 32) + (cross >> 32) + ahbh;
	return result;
}
