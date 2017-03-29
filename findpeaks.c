// Find all local maxima stronger than THRESHOLD in image A, and return a list
// of pixel locations x, y and signal levels, ordered in increasing strength.
// 

#include "__base.h"
#include "findpeaks.h"
#include "groupsort.h"
#include "strings.h"

int findpeaks(float* a, int nx, int ny, float threshold,
    float* x, float* y, float* s, int maxsources)
{
    int nsources = 0;

    float* ya = a + nx;
    for (int j = 1; j < ny - 1; j++, ya += nx)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            if (ya[i] > threshold)
            {
                float* yia = ya + i - nx;
                float b = yia[-1];

#define UPDATEB(_i) if (yia[_i] > b) b = yia[_i]
                UPDATEB(0);
                UPDATEB(1);

                yia += nx;
                UPDATEB(-1);
                UPDATEB(1);

                yia += nx;
                UPDATEB(-1);
                UPDATEB(0);
                UPDATEB(1);
#undef UPDATEB

                if (ya[i] > b)
                {
                    if (nsources < maxsources)
                    {
                        x[nsources] = (float)i + 1;
                        y[nsources] = (float)j + 1;
                        s[nsources] = ya[i];
                        nsources++;
                    }
                }
            }
        }
    }

    if (nsources == 0)
        return 0;
    if (nsources == maxsources)
        printf("FINDPEAKS: Maximum source count exceeded.\n");

    // Order in increasing strength.
    groupsort(nsources, s, x, y);
    return nsources;
}