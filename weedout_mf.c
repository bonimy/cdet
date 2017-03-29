//Weed out likely duplicate sources, i.e.those within a certain critical
//radius of a previously - seen source.The source list comes in ascending
//order of source strength and leaves in descending order.

#include "__base.h"
#include "weedout_mf.h"

int weedout_mf(float* x, float* y, float* s, int nsources, float rcrit, float pixel, BOOL* notseen,
    int* satmask, float* rsatset, int nx, int ny, int nbands, int isingleband)
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

            float* satband = rsatset + ((j * nx) + i) * nbands;
            float rcritx = isingleband == 0 ?
                max(max(max(satband[0], satband[1]), satband[2]), satband[3]) :
                satband[isingleband - 1];
            rcritx *= 2 / (3600 * pixel);

            if (satmask[(j * nx) + i] == 0 || rcritx < rcrit)
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
                float rcritxsq = rcritx * rcritx;
                int imbx = lroundf(rcritx);
                int ilo = max(i - imbx, 0);
                int ihi = min(i + imbx, nx - 1);
                int jlo = max(j - imbx, 0);
                int jhi = min(j + imbx, ny - 1);

                for (int jj = jlo; jj <= jhi; jj++)
                {
                    for (int ii = ilo; ii <= ihi; ii++)
                    {
                        if ((ii - i) * (ii - i) + (jj - j) * (jj - j) <= rcritxsq)
                            notseen[jj * nx + ii] = FALSE;
                    }
                }
            }
        }
    }

    free(mbox);

    memcpy(x, xw, m * sizeof(float));
    memcpy(y, yw, m * sizeof(float));
    memcpy(s, sw, m * sizeof(float));

    free(xw);
    free(yw);
    free(sw);

    return m;
}