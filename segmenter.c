//Segment a saturated - pixel mask into objects, each characterized by
//a location and number of saturated pixels.The criterion for membership
//is to be within a radius ISEARCHRAD of an object.

#include "__base.h"
#include "segmenter.h"

void segmenter(int* msk, int nx, int ny, int iblo, int ibhi, int jblo, int jbhi, int sz,
    int *ixsat, int *iysat, float* satradius, int* nsatpix)
{
	const int isatvalue = 100;
	const int minsat = 9;

	int irsq = sz * sz;

	int jplo = jblo + sz;
	int jphi = jbhi - sz;
	int iplo = iblo + sz;
	int iphi = ibhi - sz;

	int ipw = iphi - iplo + 1;
	int iph = jphi - jplo + 1;

	int* objmsk = (int*)malloc(ipw * iph * sizeof(int));
	memset(objmsk, 0, ipw * iph * sizeof(int));

	*ixsat = 0;
	*iysat = 0;
	*satradius = 0;

	int npixels = 0;

	int oindex = 0;
	int oindexbest = 0;
	int* src = msk + jplo * nx;
	int* dest = objmsk;

	for (int j = jplo; j <= jphi; j++, src += nx, dest += ipw)
	{
		for (int i = iplo, k = 0; i <= iphi; i++, k++)
		{
			if (src[i] == isatvalue && dest[k] == 0)
			{
				//Found a saturation source. Now we fill its neighbors.
				dest[k] = ++oindex;
				int npx = 1 + fillneighbors(msk, objmsk, nx, ny, iblo, ibhi, jblo, jbhi, sz, i, j, oindex);
				if (npixels < npx)
				{
					npixels = npx;
					oindexbest = oindex;
				}
			}
		}
	}

	if (npixels >= minsat)
	{
		src = msk + jplo * nx;
		dest = objmsk;
		int xsum = 0, ysum = 0;
		for (int j = jplo; j <= jphi; j++, src += nx, dest += ipw)
		{
			for (int i = iplo, k = 0; i <= iphi; i++, k++)
			{
				if (dest[k] == oindexbest)
				{
					xsum += i;
					ysum += j;
				}
			}
		}

		src = msk + jplo * nx;
		dest = objmsk;
		float xsat = (float)xsum / npixels;
		float ysat = (float)ysum / npixels;
		float rsum = 0;
		for (int j = jplo; j <= jphi; j++, src += nx, dest += ipw)
		{
			float ysq = j - ysat;
			ysq *= ysq;
			for (int i = iplo, k = 0; i <= iphi; i++, k++)
				if (dest[k] == oindexbest)
					rsum += (i - xsat) * (i - xsat) + ysq;
		}

		*satradius = (float)sqrt(rsum / npixels);
		*ixsat = 1 + lroundf(xsat);
		*iysat = 1 + lroundf(ysat);
		*nsatpix = npixels;
	}
}

int fillneighbors(int* msk, int* objmsk, int nx, int ny, int iblo, int ibhi, int jblo, int jbhi,
	int sz, int x, int y, int oindex)
{
	const int isatvalue = 100;
	const int minsat = 9;

	int npixels = 0;

	int irsq = sz * sz;

	int jplo = jblo + sz;
	int jphi = jbhi - sz;
	int iplo = iblo + sz;
	int iphi = ibhi - sz;

	int ipw = iphi - iplo + 1;
	int iph = jphi - jplo + 1;

	int* src = msk + (y - sz) * nx;
	int* dest = objmsk + (y - jplo - sz) * ipw + (x - iplo - sz);

	for (int j = y - sz; j <= y + sz; j++, src += nx, dest += ipw)
	{
		for (int i = x - sz, k = 0; i <= x + sz; i++, k++)
		{
			if (i < iplo || i > iphi || j < jplo || j > jphi)
				continue;
			if ((i - x) * (i - x) + (j - y) * (j - y) > irsq)
				continue;

			if (src[i] == isatvalue && dest[k] == 0)
			{
				dest[k] = oindex;
				npixels++;
				npixels += fillneighbors(msk, objmsk, nx, ny, iblo, ibhi, jblo, jbhi, sz, i, j, oindex);
			}
		}
	}

	return npixels;
}