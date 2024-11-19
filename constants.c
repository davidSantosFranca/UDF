#include "udf.h"

real k1 = 1E11;   // m^3.kmol^-1.s^-1 
real k3f = 5.9E9; // m^3.kmol^-1.s^-1
real k3b = 7.5E6; // s^-1
// mixture species I, IO3, H, I3, H2BO3, H3BO3, I2, H2O
real chargeNumber[] = {1, 1, 1, 1, 1, 0, 0, 0};