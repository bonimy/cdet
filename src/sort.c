// Sorts an array RA of length N into ascending numerical order using the
// Heapsort algorithm.N is input; RA is replaced by its sorted rearrangement.

#include "__base.h"
#include "sort.h"

void sort(int n, float* ra)
{
    float rra;

    if (n <= 1)
        return;
    int L = n / 2;
    int ir = n - 1;

loop1:
    if (L > 0)
    {
        L--;
        rra = ra[L];
    }
    else
    {
        rra = ra[ir];

        ra[ir] = ra[0];
        if (--ir == 0)
        {
            ra[0] = rra;
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
            i = j;
            j = 2 * j + 1;
        }
        else
            j = ir + 1;
        goto loop2;
    }
    ra[i] = rra;
    goto loop1;
}