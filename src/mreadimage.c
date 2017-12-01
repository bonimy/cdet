// Read a FITS image.This is the subroutine from T.Jarrett's file
// readfits.f, with minor modifications for use with MDET.

#include "__base.h"
#include "mreadimage.h"
#include "strings.h"
#include "fitsio.h"

int mreadimage(int* nx, int* ny, float* array, char* fname, double* crval1, double* crval2,
    double* cdelt1, double* cdelt2, double* crot, double* crpix1, double* crpix2)
{
    fitsfile* fptr;

    // The STATUS parameter must always be initialized.
    int status = 0;

    // Open the FITS file
    if (ffopen(&fptr, fname, READONLY, &status))
    {
        printf("MREADIMAGE: Could not read %s\n", fname);
        return 3;
    }

    // Determine the size of the image.
    long naxes[2];
    int nfound;
    ffgknj(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);

    // Check that it found both NAXIS1 and NAXIS2 keywords.
    if (nfound != 2)
    {
        printf("MREADIMAGE: Failed to read the NAXISn keywords.");
        return 4;
    }

    *nx = (int)naxes[0];
    *ny = (int)naxes[1];

    char comment[80];
    *crval1 = 0;
    if (ffgkyd(fptr, "CRVAL1", crval1, comment, &status) > 0)
    {
        *crval1 = 0;
        status = 0;
    }

    *crval2 = 0;
    if (ffgkyd(fptr, "CRVAL2", crval2, comment, &status) > 0)
    {
        *crval2 = 0;
        status = 0;
    }

    if (ffgkyd(fptr, "CRPIX1", crpix1, comment, &status) > 0)
        status = 0;
    if (ffgkyd(fptr, "CRPIX2", crpix2, comment, &status) > 0)
        status = 0;

    *cdelt1 = *cdelt2 = 0;
    if (ffgkyd(fptr, "CDELT1", cdelt1, comment, &status) > 0)
        status = 0;
    if (ffgkyd(fptr, "CDELT2", cdelt2, comment, &status) > 0)
        status = 0;

    if (ffgkyd(fptr, "CROTA2", crot, comment, &status) > 0)
    {
        *crot = 0;
        status = 0;
    }

    // CD matrix
    if (fabs(*cdelt1) <= 1e-5 && fabs(*cdelt2) <= 1e-5)
    {
        float cd1_1 = 0, cd1_2 = 0, cd2_1 = 0, cd2_2 = 0;

        ffgkye(fptr, "CD1_1", &cd1_1, comment, &status);
        status = 0;

        ffgkye(fptr, "CD1_2", &cd1_2, comment, &status);
        status = 0;

        ffgkye(fptr, "CD2_1", &cd2_1, comment, &status);
        status = 0;

        ffgkye(fptr, "CD2_2", &cd2_2, comment, &status);
        status = 0;

        float cd1a, cd1b, cd2a, cd2b;
        if (cd2_2 != 0)
        {
            float rat = cd1_2 / cd2_2;
            float angle = -atanf(rat) * 57.2957795f;

            float tdb = angle / 57.2957795f;

            cd2a = cd2_2 / cosf(tdb);
            cd2b = -cd1_2 / sinf(tdb);

            if (fabsf(cd2_2) > fabsf(cd1_2))
                *cdelt2 = cd2a;
            else
                *cdelt2 = cd2b;

            cd1a = -cd1_1 / cosf(tdb);
            cd1b = -cd2_1 / sinf(tdb);

            if (fabsf(cd1_1) > fabsf(cd2_1))
                *cdelt1 = cd1a;
            else
                *cdelt1 = cd1a;

            *crot = angle;
        }
        else if (cd1_1 != 0)
        {
            float rat = cd2_1 / cd1_1;
            float angle = -atanf(rat) * 57.2957795f;

            float tdb = angle / 57.2957795f;

            cd2a = cd2_2 / cosf(tdb);
            cd2b = -cd1_2 / sinf(tdb);

            if (fabsf(cd2_2) > fabsf(cd1_2))
                *cdelt2 = cd2a;
            else
                *cdelt2 = cd2b;

            cd1a = -cd1_1 / cosf(tdb);
            cd1b = -cd2_1 / sinf(tdb);

            if (fabsf(cd1_1) > fabsf(cd2_1))
                *cdelt1 = cd1a;
            else
                *cdelt1 = cd1a;

            *crot = angle;
        }
    }

    // Initialize variables
    int npixels = naxes[0] * naxes[1];
    int anynul;
    ffgpve(fptr, 1, 1, npixels, -999, array, &anynul, &status);

    // The FITS file must always be closed before exiting the program.
    ffclos(fptr, &status);
    return status;
}