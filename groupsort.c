// Sorts an array RA of length N into ascending numerical order using the
// Heapsort algorithm.N is input; RA is replaced by its sorted rearrangement.
// X and Y are arrays associated with RA and which get sorted along with
// RA.

#include "__base.h"
#include "groupsort.h"

void groupsort(int n, float* ra, float* x, float* y)
{
	float rra, rx, ry;

	if (n <= 1)
		return;
	int L = n / 2;
	int ir = n - 1;

loop1:
	if (L > 0)
	{
		L--;
		rra = ra[L];
		rx = x[L];
		ry = y[L];
	}
	else
	{
		rra = ra[ir];
		rx = x[ir];
		ry = y[ir];

		ra[ir] = ra[0];
		x[ir] = x[0];
		y[ir] = y[0];
		if (--ir == 0)
		{
			ra[0] = rra;
			x[0] = rx;
			y[0] = ry;
			return;
		}
	}
	int i = L;
	int j = 2 * L + 1;
loop2:
	if (j <= ir)
	{
		if (j < ir)
		{
			if (ra[j] < ra[j + 1])
				j++;
		}
		if (rra < ra[j])
		{
			ra[i] = ra[j];
			x[i] = x[j];
			y[i] = y[j];
			i = j;
			j = 2 * j + 1;
		}
		else
			j = ir + 1;
		goto loop2;
	}
	ra[i] = rra;
	x[i] = rx;
	y[i] = ry;
	goto loop1;
}