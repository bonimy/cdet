//Estimate slowly varying background of image, A, based on median - filtering
//window of specified width[pixels].Also do confusion estimation.
//
//Normally, the input parameter ISKIP should be set equal to 1.  But to make
//a quick estimate at the expense of some accuracy, set iskip at some larger
//integer, e.g. 8.

#include "__base.h"
#include "svb.h"
#include "svb_mf.h"
#include "estmode.h"
#include "gausscalc.h"
#include "medsort.h"

void svb_mf(float* a, float* unca, int nx, int ny, int width, int iskip, float* b, float* csig, BOOL multiframe)
{
    // A value indicating a blank, unused pixel.
    const float blankpix = -999;

    int hw = width / 2;

    // Center position of background
    int ix = nx / 2;
    int iy = ny / 2;

    if (iskip * 8 > width)
        iskip = 1;

    // Get the center position of the sample image
    int isampx = ix / hw;
    int isampy = iy / hw;

    // Get the size of the sample image
    int nsampx = 2 * isampx + 1;
    int nsampy = 2 * isampx + 1;
    int nsamparea = nsampx * nsampy;
    int nsampsz = nsamparea * sizeof(float);

    float* vmed = (float*)malloc(nsampsz);
    float* vcon = (float*)malloc(nsampsz);

    // Calculate medians on a coarse grid
    float* lbox = (float*)malloc(width * width * sizeof(float));
    for (int j = -isampy, n = 0; j <= isampy; j++)
    {
        int jvalue = iy + j * hw;
        for (int i = -isampx; i <= isampx; i++, n++)
        {
            int ivalue = ix + i * hw;

            float usum;
            int k = fillbox(a, unca, nx, ny, lbox, width, ivalue, jvalue, iskip, &usum);

            float q1, q2;
            float bmed = medsortdirect(lbox, k, &q1, &q2);
            float bsig = bmed - q1;
            float bmode = estmode(lbox, k, bmed, bsig);
            float uncmean = usum / k;

            vmed[n] = bmed;
            float sigpix = multiframe ? (bmode - q1) : ((q2 - q1) / 2);
            vcon[n] = (float)sqrt(max(sigpix * sigpix - uncmean * uncmean, 0));
        }
    }
    free(lbox);

    float bsum0 = medsort(vmed, nsamparea, 0, 0);
    float csum0 = medsort(vcon, nsamparea, 0, 0);

    // Interpolate using a Gaussian kernel.
    float fwhm = 1.5f * hw;
    int ik = 3 * hw / 2;
    int ksize = 2 * ik + 1;

    int ksz = ksize * ksize * sizeof(float);
    float* kernel = (float*)malloc(ksz);
    gausscalc(kernel, ksize, fwhm);

    float lowest = kernel[ik];
    for (int i = ksize * ksize; --i >= 0; )
    {
        if (kernel[i] < lowest)
            kernel[i] = 0;
        else
            kernel[i] -= lowest;
    }

    int narea = nx * ny;
    int nsz = narea * sizeof(float);
    float* bwt = (float*)malloc(nsz);
    memset(bwt, 0, nsz);

    float* vm = vmed;
    float* vc = vcon;
    float* bdest = b;
    float* cdest = csig;
    memset(b, 0, nsz);
    memset(csig, 0, nsz);
    int yskip = nx * hw;

    for (int j = -isampy; j <= isampy; j++, vm += nsampx, vc += nsampx)
    {
        int jvalue = iy + j * hw;
        for (int i = -isampx; i <= isampx; i++)
        {
            int ivalue = ix + i * hw;
            interpolate(b, csig, bwt, nx, ny, kernel, ik, vm[i + isampx], vc[i + isampx], ivalue, jvalue);
        }
    }

    free(vmed);
    free(vcon);

    for (int i = narea; --i >= 0; )
    {
        if (bwt[i] != 0)
        {
            b[i] /= bwt[i];
            csig[i] /= bwt[i];
        }
        else
        {
            b[i] = bsum0;
            csig[i] = csum0;
        }
    }

    free(bwt);
}