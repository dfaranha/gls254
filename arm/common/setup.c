#include "setup.h"

#include "utils.h"
#include "basefield.h"

void init_components() {
	utils_init();
	precomp_inv_tables();
}

void dispose_components() {
	
}
