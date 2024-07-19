#include "udf.h"
#include "functions.c"
#include "constants.c"

real chargeNumber[] = {1, 1, 1, 1, 1, 0, 0, 0};

DEFINE_VR_RATE(vol_reac_rate, c, t, r, mw, yi, rate, rr_t)
{
    // Number of species
    int N = 8;
    // mixture species I, IO3, H, I3, H2BO3, H2O, H3BO3, I2
    real mfI = yi[0];
    real mlI = mw[0];
    real mfIO3 = yi[1];
    real mlIO3 = mw[1];
    real mfH = yi[2];
    real mlH = mw[2];
    real mfI3 = yi[3];
    real mlI3 = mw[3];
    real mfH2BO3 = yi[4];
    real mlH2BO3 = mw[4];
    real mfH3BO3 = yi[5];
    real mlH3BO3 = mw[5];
    real mfI2 = yi[6];
    real mlI2 = mw[6];
    real mfH2O = yi[7];
    real mlH2O = mw[7];

    if (!strcmp(r->name, "reaction-1"))
    {
        /* Reaction 1 - k1*[H+]*[H2BO3-]*/
        *rate = k1 * mfH / mlH * mfH2BO3 / mlH2BO3;
    }
    else if (!strcmp(r->name, "reaction-2"))
    {
        /* Reaction 2 - k2*[H+]²*[I-]²*[IO3-]*/
        real molar_fraction[N];
        calculate_molar_fraction(c, t, molar_fraction, mw, N);
        real k2 = calculate_K2(molar_fraction, chargeNumber, N);
        *rate = k2 * mfH / mlH * mfH / mlH * mfI / mlI * mfI / mlI * mfIO3 / mlIO3;
    }
    else if (!strcmp(r->name, "reaction-3"))
    {
        /* Reaction 3 - k3f*[I3-] - k3b*[I2]*[I-] */
        *rate = k3f * mfI3 / mlI3 - k3b * mfI2 / mlI2 * mfI / mlI;
    }
    else
    {
        Message("Reaction not found\n");
    }

    *rr_t = *rate;
}