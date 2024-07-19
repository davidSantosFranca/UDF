#include "udf.h"

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
static real calculate_ionic_strength(real molarConcentration[], real chargeNumber[], int size)
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

/**
 * Calculates the molar fraction of species in a cell.
 *
 * @param c The cell index.
 * @param t The thread pointer.
 * @param molar_fraction An array to store the calculated molar fractions.
 * @param molecular_weight_species An array of molecular weights for each species.
 * @param N The number of species.
 */
void calculate_molar_fraction(cell_t c, Thread *t, real *molar_fraction, real *molecular_weight_species, int N)
{
    int i;
    real mass_fraction;
    real molecular_weight_mixture = 0.0;

    for (i = 0; i < N; i++)
    {
        mass_fraction = C_YI(c, t, i);
        molecular_weight_mixture += mass_fraction * molecular_weight_species[i];
    }

    for (i = 0; i < N; i++)
    {
        mass_fraction = C_YI(c, t, i);
        molar_fraction[i] = (mass_fraction * molecular_weight_mixture) / molecular_weight_species[i];
    }
}

/**
 * Calculates the value of K2 based on the given molar concentrations, charge numbers, and size.
 *
 * @param molarConcentration An array of molar concentrations.
 * @param chargeNumber An array of charge numbers.
 * @param size The size of the arrays.
 * @return The calculated value of K2.
 */
static real calculate_K2(real molarConcentration[], real chargeNumber[], int size)
{
    real ionic_strength = calculate_ionic_strength(molarConcentration, chargeNumber, size);
    if (ionic_strength < 0.166)
    {
        return pow(10, (9.28105 - 3.664 * sqrt(ionic_strength)));
    }
    else
    {
        return pow(10, (8.383 - 1.5115 * sqrt(ionic_strength) + 0.23689 * ionic_strength));
    }
}