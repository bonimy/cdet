//---------------------------------------------------------------------------
//MDET: Multiband DETection.
//
//This program performs the detection step for multiple bands using the output
//from the image coadder AWAIC.
//
//COMMAND LINE SYNTAX :
//     . / mdet[-c] - imageN <FITS file name> -sigimageN <FITS file name>
//- cmaskN <FITS file name> -backwindow <real>
//- threshold <real> -outlist <text file name> -fwhmN <real>
//- focalpix <real> -matchout <FITS file name>
//- svboutN <FITS file name> -svnprec <sigfigs>
//
//     where N = 1, 2, 3 or 4.
//
//      If the - c option is included, the flux threshold will automatically
//     be adjusted to allow for confusion.
//
//INPUT PARAMETERS :
//     imageN = name of FITS file containing input coadded
//                       image for band N(all input images are assumed
//                       to be coregistered)
//     sigimageN = name of FITS file containing the corresponding
//                       uncertainty the coadded image
//     cmaskN = name of FITS file containing the coadd mask
//     backwindow = width of median filtering window[coadd pixels]
//     threshold = detection threshold[sigmas]
//     outlist = name of output file containing a list of
//                       candidate detections
//     fwhmN = nominal PSF FWHM for band N[arcsec].Default
//                       values: 6., 6., 9., 15. for the 4 bands, resp.
//     focalpixN = focal - plane pixel size for band N[arcsec].
//                       Default values : 2.75, 2.75, 2.75, 5.5
//     matchout = name of output matched filter image;  if not
//                       specified, then do not write out
//     svbN = name of FITS file for the output slowly - varying
//                       background image in band N; if not specified,
//                       then do not write out
//      sigfigs = if >0, # of sigfigs to keep in svb pixels
//
//OUTPUTS:
// (1) List of candidate sources(text file).The format is as follows :
//
//     The first 6 records are header lines, the first of which gives the
//     number of candidates.The subsequent records(7, ... Ncandidates + 6)
//     list the following quantities :
//
//     n, RA, Dec, SNR, Rsat1, Rsat2, Rsat3, Rsat4
//
//      where :
//       n = candidate number
//       RA = RA[deg]
//       Dec = Dec[deg]
//       SNR = detection SNR of this candidate
//       RsatN = saturation radius[arcsec], where N = 1, ...4
//
// (2) Set of slowly - varying background images(FITS format).
//     These represent median - filtered versions of the input images.
//
//EXIT CODES :
//     0 = normal termination
//     1 = images misregistered(non - fatal error; proceed anyway)
//     2 = missing input parameters
//     3 = file not found
//     4 = problem with FITS keyword values
//     5 = array allocation failure
//
//DEVELOPMENT HISTORY :
//     Written by K.A.Marsh, IPAC.
//
//     Version 0.0 : IDL prototype(mdet.pro) 2008 Feb 28
//     Version 1.0 : Fortran translation(mdet.f) 2008 Apr 24
//     Version 1.1 : Blank out bright stars during background estim. 2008 Jul 9
//     Version 2.0 : Explicit declaration of real(4)
//     Version 2.5 : Add - svnprec to command line(Tim) and increase the
//                avoidance radius to 1.3*fwhm for distinguishing
//                neighboring peaks;  2009 March 19
//     Version 2.9: Raise threshold in confused regions; 2009 June 16
//     Version 3.0: For a saturated star, take candidate position as centroid
//                of the saturated region;  2009 Oct 12
//     Version 3.5: Optimize parameters for saturated star processing.
//     Version 3.5mf : Multiframe version.
//                 Saturated star bug fix 2010 Oct 4.
//                 Change WCS variables to double precision 2010 Oct 7
//     Version 4 : Bug fix to saturated star filtering.
//              4.1; Fixed(at a cost in speed) inability to get sat.radii
//                   after masking is done.
//                   Other small performance enhancements.
//     Version 5.0: Ported code to C, optimized several routines. Program
//					 completes in less than half the time now.
//-------------------------------------------------------------------------

#include "__base.h"
#include "mdet.h"
#include "appendlist.h"
#include "findpeaks.h"
#include "headpar.h"
#include "medsort.h"
#include "mparg.h"
#include "mreadimage.h"
#include "mwimage.h"
#include "pix2eq.h"
#include "segmenter.h"
#include "shuffle.h"
#include "svb.h"
#include "svb_mf.h"
#include "weedout.h"
#include "weedout_mf.h"
#include "writeimage.h"

char s0[MAX_PATH], outlist[MAX_PATH], matchout[MAX_PATH];
char imagelist[MAX_BANDS][MAX_PATH], sigimagelist[MAX_BANDS][MAX_PATH], cmasklist[MAX_BANDS][MAX_PATH];
char svblist[MAX_BANDS][MAX_PATH];
int sigfigs;
double ra, dec, x8, y8, ra0, dec0, xz, yz, pixel8, theta;
double crval1, crval2, cdelt1, cdelt2, crota, crpix1, crpix2;
double ra0c, dec0c, pixelc, thetac, x0c, y0c, rot;
double crval1m, crval2m, cdelt1m, cdelt2m, rotm, crpix1m, crpix2m;
float fwhm[MAX_BANDS], focalpix[MAX_BANDS], satradii[MAX_BANDS], annvalues[MAX_ANN];
float xcan[MAX_SOURCES], ycan[MAX_SOURCES], scan[MAX_SOURCES];
float xband[MAX_SOURCES], yband[MAX_SOURCES], sband[MAX_SOURCES];
int bandlist[MAX_BANDS], backwindow, nbands;
float threshold;
BOOL confused, multiframe, newpeak, keepblanking;
int ixbr, iybr, nx, ny;

float impixprcf(float f, int accuracy)
{
    int exp;
    float man;
    man = frexpf(f, &exp);
    return ldexpf(rintf(man*accuracy) / accuracy, exp);
}

int mdet(int argc, char** argv)
{
    outlist[0] = '\0';
    matchout[0] = '\0';
    backwindow = 0;
    threshold = 0;
    fwhm[0] = 6;
    fwhm[1] = 6;
    fwhm[2] = 9;
    fwhm[3] = 15;
    focalpix[0] = 2.75;
    focalpix[1] = 2.75;
    focalpix[2] = 2.75;
    focalpix[3] = 5.5;
    bandlist[0] = 0;
    bandlist[1] = 0;
    bandlist[2] = 0;
    bandlist[3] = 0;
    ixbr = 0;
    iybr = 0;

    printf("MDET Version 5.0; 2017 January 08\n");
	printf("C-port by Nelson Garcia under Peter Eisenhardt.\n");

    int istat = 0;

    // Read the command arguments
    mparg(argv, argc, &nbands, &backwindow, &threshold, fwhm, focalpix, bandlist, imagelist,
        cmasklist, svblist, sigimagelist, outlist, matchout, &confused, &multiframe, &sigfigs);
    if (istat > 1)
        return istat;

    // Get some FITS header information.
    istat = headpar(&nx, &ny, imagelist[0], &ra0, &dec0, &cdelt1, &pixel8, &theta, &xz, &yz);
    float pixel = (float)pixel8;
    if (istat > 1)
        return istat;
    if (pixel < 0 || cdelt1 > 0)
    {
        printf("MDET: Negative CDELT1 and positive CDELT2 required.\n");
        return 4;
    }
    int noffsets = 0;

    if (nbands >= 2)
    {
        for (int iband = 1; iband < nbands; iband++)
        {
            int naxis1, naxis2;
            istat = headpar(&naxis1, &naxis2, imagelist[iband], &crval1, &crval2, &cdelt1,
                &cdelt2, &crota, &crpix1, &crpix2);
            if (istat > 1)
                return istat;
            if (naxis1 < nx)
                nx = naxis1;
            if (naxis2 < ny)
                ny = naxis2;
        }
        if (noffsets > 0)
            printf("\n");
    }

    float* matchset = (float*)malloc(nx * ny * nbands * sizeof(float));
    float* match = (float*)malloc(nx * ny * sizeof(float));
    float* rsatset = (float*)malloc(nx * ny * nbands * sizeof(float));
    int* satmask = (int*)malloc(nx * ny * sizeof(int));
    if (matchset == 0 || match == 0 || rsatset == 0 || satmask == 0)
        return 5;

    memset(matchset, 0, nx * ny * nbands * sizeof(float));
    memset(rsatset, 0, nx * ny * nbands * sizeof(float));
    memset(match, 0, nx * ny * sizeof(float));

    // Read input images for all bands.
    for (int iband = 0; iband < nbands; iband++)
    {
        float* C = (float*)malloc(nx * ny * sizeof(float));
        float* Cm = (float*)malloc(nx * ny * sizeof(float));
        float* Cb = (float*)malloc(nx * ny * sizeof(float));
        float* U = (float*)malloc(nx * ny * sizeof(float));
        int* Msk = (int*)malloc(nx * ny * sizeof(int));
        float* b = (float*)malloc(nx * ny * sizeof(float));
        float* bb = (float*)malloc(nx * ny * sizeof(float));
        float* bsc = (float*)malloc(nx * ny * sizeof(float));
        float* csig = (float*)malloc(nx * ny * sizeof(float));
        if (C == 0 || Cm == 0 || Cb == 0 || U == 0 || Msk == 0 ||
            b == 0 || bb == 0 || bsc == 0 || csig == 0)
            return 5;

        memset(b, 0, nx * ny * sizeof(float));
        memset(bb, 0, nx * ny * sizeof(float));
        memset(csig, 0, nx * ny * sizeof(float));

        // Read the coadded image.
        int nxc, nyc;
        istat = headpar(&nxc, &nyc, imagelist[iband], &ra0c, &dec0c, &cdelt1, &pixelc, &thetac, &x0c, &y0c);
        int Lsize = nxc * nyc;
        float* Larray = (float*)malloc(Lsize * sizeof(float));
        float* Lsubarray = (float*)malloc(Lsize * sizeof(float));
        float* Array = (float*)malloc(nxc * nyc * sizeof(float));
        if (Larray == 0 || Lsubarray == 0 || Array == 0)
            return 5;
        istat = mreadimage(&nxc, &nyc, Larray, imagelist[iband], &crval1, &crval2, &cdelt1, &cdelt2, &rot, &crpix1, &crpix2);
        if (istat > 1)
            return istat;
        shuffle(Larray, nxc, nyc, Array);
        memcpy(C, Array, nx * ny * sizeof(float));
        float peakimage = BLANK_PIXEL;
        for (int i = nx * ny; --i >= 0; )
        {
            if (peakimage < C[i])
                peakimage = C[i];
        }

        // Read uncertainty image.
        int nxu, nyu;
        istat = headpar(&nxu, &nyu, sigimagelist[iband], &crval1, &crval2, &cdelt1, &cdelt2, &rot, &crpix1, &crpix2);
        if (istat > 1)
            return istat;
        if (nxu != nxc || nyu != nyc)
        {
            printf("MDET: Coadd and uncertainty images incompatible\n");
            return 4;
        }
        istat = mreadimage(&nxu, &nyu, Larray, sigimagelist[iband], &crval1, &crval2, &cdelt1, &cdelt2, &rot, &crpix1, &crpix2);
        if (istat > 1)
            return istat;
        shuffle(Larray, nxc, nyc, Array);
        memcpy(U, Array, nx * ny * sizeof(float));
        float Umax = BLANK_PIXEL;
        for (int i = nx * ny; --i >= 0; )
        {
            if (Umax < U[i])
                Umax = U[i];
        }
        for (int i = nx * ny; --i >= 0; )
        {
            if (U[i] == BLANK_PIXEL || U[i] == 0)
                U[i] = Umax;
        }

        // Replace uncertainty image with its median.
        if (multiframe)
        {
            int noblank = 0;
            for (int i = 0; i < Lsize; i++)
            {
                if (Larray[i] != BLANK_PIXEL && Larray[i] != 0)
                    Lsubarray[noblank++] = Larray[i];
            }
            float umed = medsortdirect(Lsubarray, noblank, 0, 0);
            for (int i = nx * ny; --i >= 0; )
                U[i] = umed;
        }

        if (strlen(cmasklist[iband]))
        {
            // Read coadd mask.
            int nxm, nym;
            istat = headpar(&nxm, &nym, cmasklist[iband], &crval1m, &crval2m, &cdelt1m, &cdelt2m, &rotm, &crpix1m, &crpix2m);
            if (istat > 1)
                return istat;
            if (nxm != nxc || nym != nyc)
            {
                printf("MDET: Coadd mask incompatible with image\n");
                return 4;
            }
            istat = mreadimage(&nxm, &nym, Larray, cmasklist[iband], &crval1m, &crval2m, &cdelt1m, &cdelt2m, &rotm, &crpix1m, &crpix2m);
            if (istat > 1)
                return istat;
            shuffle(Larray, nxc, nyc, Array);
            for (int i = nx * ny; --i >= 0; )
                Msk[i] = lroundf(Array[i]);
        }
        else
        {
            printf("MDET: WARNING: No coadd mask specified for band %i. Assuming all zeroes\n", iband + 1);

            memset(Msk, 0, nx * ny * sizeof(int));
        }

        // Calculate slowly-varying background and write it out.
        int width = 2 * lroundf(backwindow / 2.f) + 1;

        // First mask out the NBRIGHT stars for which SNR > BTHRESH, using a mask width
        // corresponding to half the width of the median filtering window.To detect
        // these NBRIGHT stars, use an approximate procedure in which the background
        // estimation is based on sampling every ISKIPth pixel.Limit the masked area
        // to at most a third of the image.
        memcpy(Cb, C, nx * ny * sizeof(float));
        for (int i = nx * ny; --i >= 0; )
        {
            if (Msk[i] == MAX_SAT)
                Cb[i] = BLANK_PIXEL;
        }
        int iskip = 8;
        float bthresh = 100;
        float fwhmpix = fwhm[bandlist[iband]] / (3600.f * pixel);
        int iblank = min(lroundf(25 * fwhmpix), (int)(width / 2.0f));

        if (multiframe)
            svb_mf(Cb, U, nx, ny, width, iskip, bb, csig, multiframe);
        else
            svb(Cb, U, nx, ny, width, iskip, bb, csig);

        float* peaks = (float*)malloc(nx * ny * sizeof(float));
        for (int i = nx * ny; --i >= 0; )
        {
            peaks[i] = (C[i] - bb[i]) / U[i];
        }
        int nbright = findpeaks(peaks, nx, ny, bthresh, xcan, ycan, scan, MAX_SOURCES);
        free(peaks);

        float *ibstar = (float*)malloc(nbright * sizeof(float));
        float *jbstar = (float*)malloc(nbright * sizeof(float));
        float *srbstar = (float*)malloc(nbright * sizeof(float));
        float *pbstar = (float*)malloc(nbright * sizeof(float));
        memcpy(Cm, C, nx * ny * sizeof(float));
        memcpy(Cb, C, nx * ny * sizeof(float));
        memset(bsc, 0, nx * ny * sizeof(float));
        memset(satmask, 0, nx * ny * sizeof(float));
        int ibright = 0;
        memset(srbstar, 0, nbright * sizeof(float));
        int nsatflagged = 0;
        int nsatr = 0;
        int totsat = 0;
        for (int i = nx * ny; --i >= 0; )
        {
            if (Msk[i] == MAX_SAT)
                totsat++;
        }
        int isearchrad = lroundf(focalpix[bandlist[iband]] / (3600 * pixel));

        printf("\n");
        printf("Band=%1i, NBrt=%6i, TotSat=%8.1f\n", 1 + bandlist[iband], nbright, (float)totsat);

        if (nbright != 0)
        {
            keepblanking = TRUE;
            for (int jpeak = 0; jpeak < nbright; jpeak++)
            {
                int ipeak = (nbright - 1) - jpeak;
                int ibr = lroundf(xcan[ipeak]) - 1;
                int jbr = lroundf(ycan[ipeak]) - 1;
                newpeak = FALSE;
                if (Cb[(jbr * nx) + ibr] != BLANK_PIXEL)
                    newpeak = TRUE;

                int iblo = max(ibr - iblank, 0);
                int ibhi = min(ibr + iblank, nx - 1);
                int jblo = max(jbr - iblank, 0);
                int jbhi = min(jbr + iblank, ny - 1);

                if (keepblanking)
                {
                    for (int j = jblo; j <= jbhi; j++)
                    {
                        for (int i = iblo; i <= ibhi; i++)
                            Cm[(j * nx) + i] = BLANK_PIXEL;
                    }

                    float blankpixcount = 0;
                    for (int k = nx * ny; --k >= 0; )
                    {
                        if (Cm[k] == BLANK_PIXEL)
                            blankpixcount += Cm[k];
                    }
                    blankpixcount /= BLANK_PIXEL;
                    float blankfrac = (float)blankpixcount / (nx * ny);
                    if (blankfrac > 0.3333)
                        keepblanking = FALSE;
                }

                if (newpeak)
                {
                    // blank out an area around this peak so it won't be used again
                    float blankr = 5 * fwhmpix;
                    int jblankr = lroundf(blankr);
                    float blankrsq = blankr * blankr;
                    int ilo = max(ibr - jblankr, 0);
                    int ihi = min(ibr + jblankr, nx - 1);
                    int jlo = max(jbr - jblankr, 0);
                    int jhi = min(jbr + jblankr, ny - 1);
                    int nsat = 0;
                    for (int j = jlo; j <= jhi; j++)
                    {
                        for (int i = ilo; i <= ihi; i++)
                        {
                            if ((i - ibr) * (i - ibr) + (j - jbr) * (j - jbr) <= blankrsq)
                                Cb[j * nx + i] = BLANK_PIXEL;
                            if (Msk[j * nx + i] == MAX_SAT)
                                nsat++;
                        }
                    }

                    int nsatpix = 0;
                    int ixsat, iysat;
                    float satradius;
                    if (keepblanking)
                    {
                        // get sat radius

                        // Get the regional peak
                        float cmax = -1e35f;
                        for (int j = jblo; j <= jbhi; j++)
                        {
                            for (int i = iblo; i <= ibhi; i++)
                            {
                                int index = j * nx + i;
                                if (C[index] - bb[index] > cmax)
                                {
                                    ixbr = i;
                                    iybr = j;
                                    cmax = C[index] - bb[index];
                                }
                            }
                        }

                        ibstar[ibright] = (float)ixbr + 1;
                        jbstar[ibright] = (float)iybr + 1;
                        pbstar[ibright] = cmax;
                        ibright++;

                        segmenter(Msk, nx, ny, iblo + 1, ibhi + 1, jblo + 1, jbhi + 1, isearchrad, &ixsat, &iysat, &satradius, &nsatpix);
                    }
                    else
                    {
                        // Special call on smaller area after masking is done
                        // to compute sat radius for this source.Only run if
                        // it's possible for the result to be big enough for the
                        // worth the trouble
                        if (nsat > 3.14*(fwhmpix * 0.75)*(fwhmpix * 0.75))
                        {
                            segmenter(Msk, nx, ny, iblo + 1, ibhi + 1, jblo + 1, jbhi + 1, isearchrad, &ixsat, &iysat, &satradius, &nsatpix);
                        }
                    }

                    if (nsatpix > 0)
                    {
                        nsatr++;
                        int isatradius = lroundf(satradius);
                        float rsatarcsec = satradius * (float)pixelc * 3600;
                        float dsatsq = (2 * satradius) * (2 * satradius);
                        ilo = max(ixsat - 2 * isatradius, 0);
                        ihi = min(ixsat + 2 * isatradius, nx - 1);
                        jlo = max(iysat - 2 * isatradius, 0);
                        jhi = min(iysat + 2 * isatradius, ny - 1);
                        for (int j = jblo; j <= jbhi; j++)
                        {
                            for (int i = iblo; i <= ibhi; i++)
                            {
                                int index = (j - 1) * nx + (i - 1);
                                if ((i - ixsat) * (i - ixsat) + (j - iysat) * (j - iysat) <= dsatsq)
                                {
                                    if (satmask[index] == 0)
                                    {
                                        if (Msk[index] == MAX_SAT)
                                            nsatflagged++;
                                        satmask[index] = 1;
                                        rsatset[iband * nx * ny + index] = rsatarcsec;
                                    }
                                }
                            }
                        }

                        if (keepblanking)
                        {
                            ibstar[ibright - 1] = (float)ixsat;
                            jbstar[ibright - 1] = (float)iysat;
                            srbstar[ibright - 1] = satradius;
                        }
                    }
                }
            }
        }

        printf("Number of non-zero saturation radii = %i\n", nsatr);

        // Also blank out all areas flagged with saturation in prep for svb computation ! TPC
        // ... or not
        for (int i = nx * ny; --i >= 0; )
        {
            if (Msk[i] == MAX_SAT)
                Cm[i] = BLANK_PIXEL;
        }

        // Now median filter the remaining pixels.
        printf("Computing SVB for band %i\n", bandlist[iband] + 1);
        if (multiframe)
            svb_mf(Cm, U, nx, ny, width, 1, b, csig, multiframe);
        else
            svb(Cm, U, nx, ny, width, 1, b, csig);

        // Increase the local confusion error in the vicinity of each saturated star.
        if (ibright > 0)
        {
            float annwid = fwhmpix / 4;
            float rin = 3 * fwhmpix;
            int nannuli = (int)((iblank - rin) / annwid);

            for (int ii = 0; ii < ibright; ii++)
            {
                int ixbr = (int)ibstar[ii] - 1;
                int iybr = (int)jbstar[ii] - 1;
                if (satmask[iybr * nx + ixbr] == 1)
                {
                    int iblo = max(ixbr - iblank, 0);
                    int ibhi = min(ixbr + iblank, nx - 1);
                    int jblo = max(iybr - iblank, 0);
                    int jbhi = min(iybr + iblank, ny - 1);
                    for (int j = jblo; j <= jbhi; j++)
                    {
                        for (int i = iblo; i <= ibhi; i++)
                        {
                            float rval = (float)sqrt((i - ixbr) * (i - ixbr) + (j - iybr) * (j - iybr));
                            if (rval < rin)
                                bsc[j * nx + i] = pbstar[ii];
                        }
                    }
                    for (int ian = 0; ian < nannuli; ian++)
                    {
                        float ran = rin + (ian)* annwid;
                        int nvalues = 0;
                        for (int j = jblo; j <= jbhi; j++)
                        {
                            for (int i = iblo; i <= ibhi; i++)
                            {
                                int index = j * nx + i;
                                float rval = (float)sqrt((i - ixbr) * (i - ixbr) + (j - iybr) * (j - iybr));
                                if (rval > ran && rval <= ran + annwid && nvalues < MAX_ANN)
                                    annvalues[nvalues++] = C[index] - b[index];
                            }
                        }
                        if (nvalues > 4)
                        {
                            float q1, q2;
                            float amed = medsort(annvalues, nvalues, &q1, &q2);
                            float annsig = (q2 - q1) / 2;
                            for (int j = jblo; j <= jbhi; j++)
                            {
                                for (int i = iblo; i <= ibhi; i++)
                                {
                                    int index = j * nx + i;
                                    float rval = (float)sqrt((i - ixbr) * (i - ixbr) + (j - iybr) * (j - iybr));
                                    if (rval > ran && rval <= ran + annwid)
                                        bsc[index] = (float)sqrt(amed * amed + annsig * annsig);
                                }
                            }
                        }
                        else
                        {
                            for (int j = jblo; j <= jbhi; j++)
                            {
                                for (int i = iblo; i <= ibhi; i++)
                                    bsc[j * nx + i] = pbstar[ii];
                            }
                        }
                    }
                    bsc[iybr * nx + ixbr] = 0;
                    C[iybr * nx + ixbr] += peakimage * srbstar[ii] * srbstar[ii];
                    C[iybr * nx + ixbr] *= 100.f;
                }
            }
        }

        free(ibstar);
        free(jbstar);
        free(srbstar);
        free(pbstar);

        // Add this term to the sum.
        if (!confused)
            memset(csig, 0, nx * ny * sizeof(float));
        for (int i = nx * ny; --i >= 0; )
        {
            if (C[i] != BLANK_PIXEL)
                C[i] -= b[i];
            if (C[i] < 0)
                C[i] = 0;

            if (bsc[i] > csig[i])
                csig[i] = bsc[i];

            C[i] = C[i] * C[i] / (U[i] * U[i] + csig[i] * csig[i]);
            match[i] += C[i];

            matchset[(nx * ny * iband) + i] = C[i];
        }

        // If requested, reduce pixel precision
        if (svblist[iband][0] != '\0')
        {
            if (sigfigs > 0)
            {
                int iscale = (int)pow(2, sigfigs / log10(2) + 1);
                for (int i = nx * ny; --i >= 0; )
                    b[i] = impixprcf(b[i], iscale);
            }
            writeimage(b, nx, ny, imagelist[0], svblist[iband]);
        }
        free(Larray);
        free(Lsubarray);
        free(Array);
        free(C);
        free(Cm);
        free(Cb);
        free(U);
        free(Msk);
        free(b);
        free(bb);
        free(bsc);
        free(csig);
    }

    // Calculate multiband matched filter image and its set of single - band
    // counterparts.
    for (int i = nx * ny; --i >= 0; )
    {
        match[i] = sqrtf(match[i]);
    }
    for (int i = nx * ny * nbands; --i >= 0; )
    {
        matchset[i] = sqrtf(matchset[i]);
    }

    // Write out the matched filter image if necessary.
    if (matchout[0] != '\0')
        writeimage(match, nx, ny, imagelist[0], matchout);

    // Find candidate sources.
    int nsources = findpeaks(match, nx, ny, threshold, xcan, ycan, scan, MAX_SOURCES);

    // Weed out duplicates, i.e.those within a certain critical radius of each
    // other.
    if (nsources > 0)
    {
        BOOL* notseen = (BOOL*)malloc(nx * ny * sizeof(BOOL));
        if (notseen == 0)
            return 5;
        memset(notseen, TRUE, nx * ny * sizeof(BOOL));
        float rcrit = min(min(min(fwhm[0], fwhm[1]), fwhm[2]), fwhm[3]);
        rcrit *= 1.3f / (3600.f * pixel);
        memset(satmask, 0, nx * ny * sizeof(int));
        for (int iband = 0; iband < nbands; iband++)
        {
            for (int i = nx * ny; --i >= 0; )
                satmask[i] += lroundf(100 * rsatset[(nx * ny * iband) + i]);
        }
        if (multiframe)
            nsources = weedout_mf(xcan, ycan, scan, nsources, rcrit, pixel, notseen, satmask, rsatset, nx, ny, nbands, 0);
        else
            nsources = weedout(xcan, ycan, scan, nsources, rcrit, notseen, satmask, nx, ny);

        // Append the list with the results of single-band detections.
        if (nbands > 1)
        {
            for (int iband = 0; iband < nbands; iband++)
            {
                memcpy(match, matchset + nx * ny * iband, nx * ny * sizeof(float));
                int newsources = findpeaks(match, nx, ny, threshold, xband, yband, sband, MAX_SOURCES);
                if (newsources > 0)
                {
                    // First, regenerate the "notseen" array.
                    memset(notseen, TRUE, nx * ny * sizeof(BOOL));
                    float rcrit = fwhm[bandlist[iband]];
                    rcrit *= 1.3f / (3600.f * pixel);
                    int imb = lroundf(rcrit);
                    int nmb = 2 * imb + 1;
                    BOOL* mbox = (BOOL*)malloc(nmb * nmb * sizeof(BOOL));
                    memset(mbox, TRUE, nmb * nmb * sizeof(BOOL));
                    float rcritsq = rcrit * rcrit;
                    for (int j = 0; j < nmb; j++)
                    {
                        for (int i = 0; i < nmb; i++)
                        {
                            if ((i - imb) * (i - imb) + (j - imb) * (j - imb) <= rcritsq)
                                mbox[j * nmb + i] = FALSE;
                        }
                    }
                    for (int n = 0; n < nsources; n++)
                    {
                        int i = lroundf(xcan[n]) - 1;
                        int j = lroundf(ycan[n]) - 1;
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
                    free(mbox);

                    // Now get rid of sources from this band which have been seen before.
                    for (int i = nx * ny; --i >= 0; )
                    {
                        satmask[i] = lroundf(100 * rsatset[(nx * ny * iband) + i]);
                    }
                    if (multiframe)
                        newsources = weedout_mf(xband, yband, sband, newsources, rcrit, pixel, notseen, satmask, rsatset, nx, ny, nbands, iband);
                    else
                        newsources = weedout(xband, yband, sband, newsources, rcrit, notseen, satmask, nx, ny);

                    // Finally, append the unique sources from this band to the master list.
                    nsources = appendlist(xcan, ycan, scan, nsources, xband, yband, sband, newsources, MAX_SOURCES);
                }
            }
        }

        free(notseen);
        printf("There were %6i candidated above SNR =%6.2f\n", nsources, threshold);
    }
    else
        printf("No candidate sources found.\n");

    FILE *f = fopen(outlist, "w");

    printf("Writing out %s\n", outlist);
    fprintf(f, "\\Nsrcs =%10i\n", nsources);
	fprintf(f, "| Src  |    RA    |    Dec   |   SNR   | Rsat1 | Rsat2 | Rsat3 | Rsat4 |\n");
	fprintf(f, "|  i   |    r     |     r    |    r    |   r   |   r   |   r   |   r   |\n");
	fprintf(f, "|      |   deg    |    deg   |         | arcsec| arcsec| arcsec| arcsec|\n");

	for (int n = 0; n < nsources; n++)
	{
		x8 = (double)xcan[n];
		y8 = (double)ycan[n];
		pix2eq(x8, y8, ra0, dec0, xz, yz, pixel8, theta, &ra, &dec);
		memset(satradii, 0, MAX_BANDS * sizeof(float));
		for (int iband = 0; iband < nbands; iband++)
		{
			float rsatband = rsatset[(iband * nx * ny) + (lroundf(ycan[n] - 1)) * nx + (lroundf(xcan[n] - 1))];
			if (rsatband > 0)
				satradii[iband] = rsatband + focalpix[iband];
		}
		while (ra < 0)
			ra += 360.;
		while (ra >= 360.)
			ra -= 360.;
		fprintf(f, "%7i%11.5f%11.5f%10.3E%8.2f%8.2f%8.2f%8.2f\n", n + 1, ra, dec, scan[n],
			satradii[0], satradii[1], satradii[2], satradii[3]);
	}

    fclose(f);
    free(match);
    free(matchset);
    free(rsatset);
    free(satmask);
    return istat;
}