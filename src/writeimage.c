// Write out a FITS file using the header from REFIMAGE.

#include "__base.h"
#include "writeimage.h"
#include "mwimage.h"

void writeimage(float* array, int nx, int ny, char* refimage, char* outfile)
{
    int lsize = nx * ny;
    float* larray = (float*)malloc(lsize * sizeof(float));
    memcpy(larray, array, lsize * sizeof(float));

    printf("Writing out %s\n", outfile);
    mwimage(nx, ny, larray, refimage, outfile);
    free(larray);
}