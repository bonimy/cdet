// Calculate a Gaussian image on a square array of NCELLS x NCELLS.

#include "__base.h"
#include "gausscalc.h"

void gausscalc(float* gauss, int ncells, float fwhm)
{
	float coef = -4 * (float)log(2) / (fwhm * fwhm);
	float* sq = (float*)malloc(ncells * sizeof(float));
	int icent = ncells / 2;

	for (int i = ncells, x = -icent; --i >= 0; x++)
		sq[i] = (float)x * x;

    float* line = gauss;
	for (int j = ncells; --j >= 0; line += ncells)
	{
		float ysq = sq[j];
		for (int i = ncells; --i >= 0; )
		{
			float arg = coef * (sq[i] + ysq);
			if (arg > -50)
				line[i] = (float)exp(arg);
		}
	}

	free(sq);
}