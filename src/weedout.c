//Weed out likely duplicate sources, i.e.those within a certain critical
//radius of a previously - seen source.The source list comes in ascending
//order of source strength and leaves in descending order.

#include "__base.h"
#include "weedout.h"

int weedout(float* x, float* y, float* s, int nsources, float rcrit, BOOL* notseen, int* satmask, int nx, int ny)
{
    int imb = lroundf(rcrit);
    int nmb = 2 * imb + 1;
    BOOL* mbox = (BOOL*)malloc(nmb * nmb * sizeof(BOOL));
    float rcritsq = rcrit * rcrit;

    BOOL* mline = mbox;
    for (int j = 0; j < nmb; j++, mline += nmb)
    {
        int ysq = j - imb;
        ysq *= ysq;

        for (int i = 0, x = -imb; i < nmb; i++, x++)
            mline[i] = rcritsq < x * x + ysq ? -1 : 0;
    }

    const float rmult = 5;
    int imbx = lroundf(rmult * rcrit);
    int nmbx = 2 * imbx + 1;
    BOOL* mbox_sat = (BOOL*)malloc(nmbx * nmbx * sizeof(BOOL));
    float rcritsqx = (rmult * rcrit) * (rmult * rcrit);

    mline = mbox_sat;
    for (int j = 0; j < nmbx; j++, mline += nmbx)
    {
        int ysq = j - imbx;
        ysq *= ysq;

        for (int i = 0, x = -imbx; i < nmbx; i++, x++)
            mline[i] = rcritsqx < x * x + ysq ? -1 : 0;
    }

    float* xw = (float*)malloc(nsources * sizeof(float));
    float* yw = (float*)malloc(nsources * sizeof(float));
    float* sw = (float*)malloc(nsources * sizeof(float));
    
    int m = 0;
    for (int nn = 0; nn < nsources; nn++)
    {
        int n = (nsources - 1) - nn;
        int i = lroundf(x[n]) - 1;
        int j = lroundf(y[n]) - 1;

        if (notseen[(j * nx) + i])
        {
            xw[m] = x[n];
            yw[m] = y[n];
            sw[m] = s[n];
            m++;

            if (satmask[(j * nx) + i] == 0)
            {
                int ilo = max(i - imb, 0);
                int ihi = min(i + imb, nx - 1);
                int jlo = max(j - imb, 0);
                int jhi = min(j + imb, ny - 1);

                for (int jj = jlo; jj <= jhi; jj++)
                {
                    for (int ii = ilo; ii <= ihi; ii++)
                    {
                        if (!mbox[(imb - j + jj) * nmb + (imb - i + ii)])
                            notseen[jj * nx + ii] = FALSE;
                    }
                }
            }
            else
            {
                int ilo = max(i - imbx, 0);
                int ihi = min(i + imbx, nx - 1);
                int jlo = max(j - imbx, 0);
                int jhi = min(j + imbx, ny - 1);

                for (int jj = jlo; jj <= jhi; jj++)
                {
                    for (int ii = ilo; ii <= ihi; ii++)
                    {
                        if (!mbox_sat[(imbx - j + jj) * nmbx + (imbx - i + ii)])
                            notseen[jj * nx + ii] = FALSE;
                    }
                }
            }
        }
    }

    free(mbox);
    free(mbox_sat);

    memcpy(x, xw, m * sizeof(float));
    memcpy(y, yw, m * sizeof(float));
    memcpy(s, sw, m * sizeof(float));

    free(xw);
    free(yw);
    free(sw);

    return m;
}