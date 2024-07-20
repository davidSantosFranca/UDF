#include "udf.h"

real calculate_ionic_strength(real molarConcentration[], real chargeNumber[], int size);
void calculate_molar_concentration(cell_t c, Thread *t, real *molar_concentration, real *molecular_weight_species, int N);
real calculate_K2(real molarConcentration[], real chargeNumber[], int size);