// Read parameter values from FITS header.

#include "__base.h"
#include "headpar.h"
#include "strings.h"
#include "fitsio.h"

int headpar(int* nx, int* ny, char* fname, double* crval1, double* crval2, double* cdelt1, double* cdelt2,
    double* crot, double* crpix1, double* crpix2)
{
    fitsfile* fptr;

    // The STATUS parameter must always be initialized.
    int status = 0;

    // Open the FITS file
    if (ffopen(&fptr, fname, READONLY, &status))
    {
        printf("HEADPAR: Could not read %s\n", fname);
        return 3;
    }

    // Determine the size of the image.
    long naxes[2];
    int nfound;
    ffgknj(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);

    // Check that it found both NAXIS1 and NAXIS2 keywords.
    if (nfound != 2)
    {
        printf("HEADPAR: Failed to read the NAXISn keywords.\n");
        return 4;
    }

    *nx = (int)naxes[0];
    *ny = (int)naxes[1];

    char comment[FLEN_COMMENT];
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

    // The FITS file must always be closed before exiting the program.
    ffclos(fptr, &status);
    return status;
}