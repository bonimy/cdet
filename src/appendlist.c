// Append the multiband detection list with the results of single-band
// detections.

#include "__base.h"
#include "appendlist.h"
#include "groupsort.h"
#include "strings.h"

int appendlist(float* x, float* y, float* s, int ncurrent,
	float* xband, float* yband, float* sband, int nband, int mx)
{
	if (ncurrent + nband > mx)
	{
		printf("APPENDLIST: Maximum source count exceeded.\n");
		nband = mx - ncurrent;
	}

    int nsz = nband * sizeof(float);
	memcpy(x + ncurrent, xband, nsz);
	memcpy(y + ncurrent, yband, nsz);
	memcpy(s + ncurrent, sband, nsz);

	ncurrent += nband;
	groupsort(ncurrent, s, x, y);

    int msz = mx * sizeof(float);
	memcpy(xband, x, msz);
	memcpy(yband, y, msz);
	memcpy(sband, s, msz);

	for (int i = ncurrent, j = 0; --i >= 0; j++)
	{
		x[j] = xband[i];
		y[j] = yband[i];
		s[j] = sband[i];
	}

	return ncurrent;
}