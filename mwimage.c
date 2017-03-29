// Create a FITS primary array containing a 2 - D image.
//
// This is T.Jarrett's subroutine wimage.f, with slight modifications 
// for use with MDET.

#include "__base.h"
#include "mwimage.h"
#include "wimage.h"
#include "fitsio.h"

void mwimage(int nx, int ny, float* larray, char* fin, char* fout)
{
    // The STATUS parameter must be initialized before using FITSIO.A
    // positive value of STATUS is returned whenever a serious error occurs.
    // FITSIO uses an `inherited status' convention, which means that if a
    // subroutine is called with a positive input value of STATUS, then the
    // subroutine will exit immediately, preserving the status value.For
    // simplicity, this program only checks the status value at the end of
    // the program, but it is usually better practice to check the status
    // value more frequently.
    int status = 0;
    deletefile(fout, &status, 0, 1);

    fitsfile* inptr;
    fitsfile* outptr;

    // The input FITS file is opened with READONLY access, and the output
    // FITS file is opened with WRITE access.
    ffopen(&inptr, fin, READONLY, &status);
    status = 0;

    // Create the new empty FITS file.The blocksize parameter is a
    // historical artifact and the value is ignored by FITSIO.
    ffinit(&outptr, fout, &status);

    // This do - loop of calls to FTGREC and FTPREC copies all the keywords from
    // the input to the output FITS file.Notice that the specified number
    // of rows in the output table, as given by the NAXIS2 keyword, will be
    // incorrect.This value will be modified later after it is known how many
    // rows will be in the table, so it does not matter how many rows are specified
    // initially.

    // Find the number of keywords in the input table header.
    int nkeys, nspace;
    ffghsp(inptr, &nkeys, &nspace, &status);

    char record[80+1];
    for (int i = 0; i < nkeys; )
    {
        ffgrec(inptr, ++i, record, &status);
        ffprec(outptr, record, &status);
    }

    // Initialize parameters about the FITS image.
    // BITPIX = 16 means that the image pixels will consist of 16 - bit
    // integers.The size of the image is given by the NAXES values.
    // The EXTEND = TRUE parameter indicates that the FITS file
    // may contain extensions following the primary array.

    // Write the array to the FITS file.
    // The last letter of the subroutine name defines the datatype of the
    // array argument; in this case the 'J' indicates that the array has an
    // integer * 4 datatype. ('I' = I * 2, 'E' = Real * 4, 'D' = Real * 8).
    // The 2D array is treated as a single 1 - D array with NAXIS1 * NAXIS2
    // total number of pixels.GROUP is seldom used parameter that should
    // almost always be set = 1.
    ffppre(outptr, 1, 1, nx * ny, larray, &status);

    // Write another optional keyword to the header
    // The keyword record will look like this in the FITS file :
    
    /*  ffdkey(outptr, "NAXIS1", &status);
        status = 0;
        ffpkyj(outptr, "NAXIS1", nx, "array size", &status);
        status = 0;
        ffdkey(outptr, "NAXIS2", &status);
        status = 0;
        ffpkyj(outptr, "NAXIS2", ny, "array size", &status);
        
        if (status != 0)
            printf("write status %i", status);
    */

    // The FITS file must always be closed before exiting the program.
    ffclos(inptr, &status);
    ffclos(outptr, &status);
}