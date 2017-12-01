// Estimate mode of the distribution of N values of an array A.

#include "__base.h"
#include "estmode.h"

#define FRSIG 25                // precision, expressed as a fraction of sigma
#define MAXSIGS 100             // maximum number of sigmas to search on
                                // either side of median
#define MAXDEVS MAXSIGS*FRSIG
#define NBINS 2*MAXDEVS

float estmode(float* a, int n, float amed, float sigma)
{
	if (n <= 0)
		return *a;

    float binsize = sigma / FRSIG;

    int h[NBINS];
	memset(h, 0, NBINS * sizeof(int));

	for (int i = n; --i >= 0; )
	{
		int ibin = MAXDEVS + lroundf((a[i] - amed) / binsize);
		if (ibin >= 0 && ibin < NBINS)
			h[ibin]++;
	}

	int ibinmode = 0;		// JWF B60519
	int npeak = -(1 << 30);
	for (int ibin = NBINS; --ibin >= 0; )
	{
		if (h[ibin] >= npeak)
		{
			npeak = h[ibin];
			ibinmode = ibin;
		}
	}

	if (npeak < 6)
		return amed;

	return amed + (ibinmode - MAXDEVS) * binsize;
}

#undef FRSIG
#undef MAXSIGS

#undef MAXDEVS
#undef NBINS