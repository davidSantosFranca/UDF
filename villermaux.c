#include "udf.h"
#include "constants.h"
#include "functions.h"

// Number of species
int n_components = 8;

DEFINE_VR_RATE(vol_reac_rate, c, t, r, mw, yi, rate, rr_t)
{
    // mixture species I, IO3, H, I3, H2BO3, H3BO3, I2, H2O
    real mfI = yi[0];
    real mwI = mw[0];     // kg/kmol
    real mfIO3 = yi[1]; 
    real mwIO3 = mw[1];   // kg/kmol
    real mfH = yi[2];
    real mwH = mw[2];     // kg/kmol
    real mfI3 = yi[3];
    real mwI3 = mw[3];    // kg/kmol
    real mfH2BO3 = yi[4];
    real mwH2BO3 = mw[4]; // kg/kmol
    real mfH3BO3 = yi[5];
    real mwH3BO3 = mw[5]; // kg/kmol
    real mfI2 = yi[6];
    real mwI2 = mw[6];    // kg/kmol
    real mfH2O = yi[7];
    real mwH2O = mw[7];   // kg/kmol

    real density = C_R(c, t); // kg/m^3

    real mcI = mfI * density / mwI;             // kmol/m^3
    real mcIO3 = mfIO3 * density / mwIO3;       // kmol/m^3
    real mcH = mfH * density / mwH;             // kmol/m^3
    real mcI3 = mfI3 * density / mwI3;          // kmol/m^3
    real mcH2BO3 = mfH2BO3 * density / mwH2BO3; // kmol/m^3
    real mcH3BO3 = mfH3BO3 * density / mwH3BO3; // kmol/m^3
    real mcI2 = mfI2 * density / mwI2;          // kmol/m^3
    real mcH2O = mfH2O * density / mwH2O;       // kmol/m^3

    if (!strcmp(r->name, "reaction-1"))
    {
        /* Reaction 1 - k1*[H+]*[H2BO3-]*/
        *rate = k1 * mcH * mcH2BO3;
    }
    else if (!strcmp(r->name, "reaction-2"))
    {
        /* Reaction 2 - k2*[H+]²*[I-]²*[IO3-]*/
        real molar_concentration[n_components]; // mol*^L-1
        calculate_molar_concentration(c, t, molar_concentration, mw, n_components);
        real k2 = calculate_K2(molar_concentration, chargeNumber, n_components); // mol*L^-1*s^-1 => kmol*m^-3*s^-1
        *rate = k2 * pow(mcH, 2) * pow(mcI, 2) * mcIO3;
    }
    else if (!strcmp(r->name, "reaction-3"))
    {
        /* Reaction 3 - k3b*[I2]*[I-] */
        *rate = k3f * mcI2 * mcI - k3b * mcI3;
    }
    else
    {
        Message("Reaction not found\n");
    }

    *rr_t = *rate;
}