#include "udf.h"
#include "functions.h"

/**
 * Function to calculate the ionic strength.
 *
 * This function calculates the ionic strength based on the molar concentrations and charge numbers of each species.
 * If the charge numbers are not provided, the function assumes all species have a charge of 1.
 *
 * @param molarConcentration Array of molar concentrations of each species.
 * @param chargeNumber Array of charge numbers of each species. If not provided (NULL), the function assumes all species have a charge of 1.
 * @param size The size of the arrays molarConcentration and chargeNumber.
 * @return The calculated ionic strength.
 */
real calculate_ionic_strength(real molarConcentration[], real chargeNumber[], int size)
{
    real sum = 0.0;
    if (chargeNumber == NULL)
    {
        for (int i = 0; i < size; i++)
        {
            sum += molarConcentration[i];
        }
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            sum += molarConcentration[i] * chargeNumber[i] * chargeNumber[i];
        }
    }
    return sum * 0.5;
}

void calculate_molar_concentration(cell_t c, Thread *t, real *molar_concentration, real *molecular_weight_species, int N)
{
    real density = C_R(c, t) * 1E-3; // kg/m^3 => kg/L
    for (int i = 0; i < N; i++)
    {
        molar_concentration[i] = C_YI(c, t, i) * density / molecular_weight_species[i];
    }
}

/**
 * Calculates the value of K2 based on the given molar concentrations, charge numbers, and size.
 *
 * @param molarConcentration An array of molar concentrations[mol/L].
 * @param chargeNumber An array of charge numbers.
 * @param size The size of the arrays.
 * @return The calculated value of K2[mol*L^-1*s^-1].
 */
real calculate_K2(real molarConcentration[], real chargeNumber[], int size)
{
    real sum = 0.0;
    if (chargeNumber == NULL)
    {
        for (int i = 0; i < size; i++)
        {
            sum += molarConcentration[i];
        }
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            sum += molarConcentration[i] * pow(chargeNumber[i], 2);
        }
    }
    real ionic_strength = sum * 0.5;
    if (ionic_strength < 0.166)
    {
        return pow(10, (9.28105 - 3.664 * sqrt(ionic_strength)));
    }
    else
    {
        return pow(10, (8.383 - 1.5115 * sqrt(ionic_strength) + 0.23689 * ionic_strength));
    }
}