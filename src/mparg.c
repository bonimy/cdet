#include "__base.h"
#include "mparg.h"

int findkey(char** argv, int argc, char* key)
{
    for (int i = argc; --i >= 0; )
    {
        if (strcmp(argv[i], key) == 0)
            return i;
    }

    return -1;
}

int mparg(char** argv, int argc, int* nbands, int* backwindow, float* threshold,
    float fwhm[MAX_BANDS], float focalpix[MAX_BANDS], BOOL bandlist[MAX_BANDS],
    char imagelist[MAX_BANDS][MAX_PATH], char cmasklist[MAX_BANDS][MAX_PATH],
    char svblist[MAX_BANDS][MAX_PATH], char sigimagelist[MAX_BANDS][MAX_PATH],
    char outlist[MAX_PATH], char matchout[MAX_PATH], BOOL* confused, BOOL* multiframe, int* sigfigs)
{
    int kimage1 = findkey(argv, argc, "-image1"),
        kimage2 = findkey(argv, argc, "-image2"),
        kimage3 = findkey(argv, argc, "-image3"),
        kimage4 = findkey(argv, argc, "-image4"),

        ksigimage1 = findkey(argv, argc, "-sigimage1"),
        ksigimage2 = findkey(argv, argc, "-sigimage2"),
        ksigimage3 = findkey(argv, argc, "-sigimage3"),
        ksigimage4 = findkey(argv, argc, "-sigimage4"),

        kcmask1 = findkey(argv, argc, "-cmask1"),
        kcmask2 = findkey(argv, argc, "-cmask2"),
        kcmask3 = findkey(argv, argc, "-cmask3"),
        kcmask4 = findkey(argv, argc, "-cmask4"),

        kbackwindow = findkey(argv, argc, "-backwindow"),
        kthreshold = findkey(argv, argc, "-threshold"),

        kfwhm1 = findkey(argv, argc, "-fwhm1"),
        kfwhm2 = findkey(argv, argc, "-fwhm2"),
        kfwhm3 = findkey(argv, argc, "-fwhm3"),
        kfwhm4 = findkey(argv, argc, "-fwhm4"),

        kfocalpix1 = findkey(argv, argc, "-focalpix1"),
        kfocalpix2 = findkey(argv, argc, "-focalpix2"),
        kfocalpix3 = findkey(argv, argc, "-focalpix3"),
        kfocalpix4 = findkey(argv, argc, "-focalpix4"),

        koutlist = findkey(argv, argc, "-outlist"),
        kmatchout = findkey(argv, argc, "-matchout"),

        ksvb1 = findkey(argv, argc, "-svb1"),
        ksvb2 = findkey(argv, argc, "-svb2"),
        ksvb3 = findkey(argv, argc, "-svb3"),
        ksvb4 = findkey(argv, argc, "-svb4"),

        ksvbprec = findkey(argv, argc, "-svbprec"),
        kc = findkey(argv, argc, "-c"),
        km = findkey(argv, argc, "-m");

#define GETKEY(key) (key == -1 ? "" : argv[key+1])

    strcpy(imagelist[0], GETKEY(kimage1));
    strcpy(imagelist[1], GETKEY(kimage2));
    strcpy(imagelist[2], GETKEY(kimage3));
    strcpy(imagelist[3], GETKEY(kimage4));

    strcpy(sigimagelist[0], GETKEY(ksigimage1));
    strcpy(sigimagelist[1], GETKEY(ksigimage2));
    strcpy(sigimagelist[2], GETKEY(ksigimage3));
    strcpy(sigimagelist[3], GETKEY(ksigimage4));

    strcpy(cmasklist[0], GETKEY(kcmask1));
    strcpy(cmasklist[1], GETKEY(kcmask2));
    strcpy(cmasklist[2], GETKEY(kcmask3));
    strcpy(cmasklist[3], GETKEY(kcmask4));
    
    if (kbackwindow != -1)
        *backwindow = (int)strtof(GETKEY(kbackwindow), 0);
    if (kthreshold != -1)
        *threshold = strtof(GETKEY(kthreshold), 0);

    if (kfwhm1 != -1)
        fwhm[0] = strtof(GETKEY(kfwhm1), 0);
    if (kfwhm2 != -1)
        fwhm[1] = strtof(GETKEY(kfwhm2), 0);
    if (kfwhm3 != -1)
        fwhm[2] = strtof(GETKEY(kfwhm3), 0);
    if (kfwhm4 != -1)
        fwhm[3] = strtof(GETKEY(kfwhm4), 0);

    if (kfocalpix1 != -1)
        focalpix[0] = strtof(GETKEY(kfocalpix1), 0);
    if (kfocalpix2 != -1)
        focalpix[1] = strtof(GETKEY(kfocalpix2), 0);
    if (kfocalpix3 != -1)
        focalpix[2] = strtof(GETKEY(kfocalpix3), 0);
    if (kfocalpix4 != -1)
        focalpix[3] = strtof(GETKEY(kfocalpix4), 0);

    strcpy(outlist, GETKEY(koutlist));
    strcpy(matchout, GETKEY(kmatchout));

    strcpy(svblist[0], GETKEY(ksvb1));
    strcpy(svblist[1], GETKEY(ksvb2));
    strcpy(svblist[2], GETKEY(ksvb3));
    strcpy(svblist[3], GETKEY(ksvb4));

    if (ksvbprec != -1)
        *sigfigs = lroundf(strtof(GETKEY(ksvbprec), 0));
    *confused = kc != -1;
    *multiframe = km != -1;

#undef GETKEY

    printf("\nINPUT PARAMETERS:\n");
    printf("image1 = %s\n", imagelist[0]);
    printf("image2 = %s\n", imagelist[1]);
    printf("image3 = %s\n", imagelist[2]);
    printf("image4 = %s\n", imagelist[3]);
    printf("sigimage1 = %s\n", sigimagelist[0]);
    printf("sigimage2 = %s\n", sigimagelist[1]);
    printf("sigimage3 = %s\n", sigimagelist[2]);
    printf("sigimage4 = %s\n", sigimagelist[3]);
    printf("backwindow = %6.1f coadd pixels\n", (float)*backwindow);
    printf("threshold = %6.2f sigmas\n", *threshold);
    printf("fwhm1-4 = %6.2f%6.2f%6.2f%6.2f arcsec\n", fwhm[0], fwhm[1], fwhm[2], fwhm[3]);
    printf("focalpix1-4 = %6.2f%6.2f%6.2f%6.2f arcsec\n", focalpix[0], focalpix[1], focalpix[2], focalpix[3]);
    printf("outlist = %s\n", outlist);
    printf("matchout = %s\n", matchout);
    printf("svb1 = %s\n", svblist[0]);
    printf("svb2 = %s\n", svblist[1]);
    printf("svb3 = %s\n", svblist[2]);
    printf("svb4 = %s\n", svblist[3]);
    printf("SVB sigfigs = %2i\n", *sigfigs);
    if (confused)
        printf("Flux threshold will be adjusted for confusion\n");
    if (multiframe)
        printf("USING MULTIFRAME VERSION\n");
    printf("\n");

    int k[] = {kimage1, kimage2, kimage3, kimage4};
    for (int i = 0; i < 4; i++)
    {
        if (k[i] != -1)
        {
            if (i != *nbands)
            {
                strcpy(imagelist[*nbands], imagelist[i]);
                strcpy(sigimagelist[*nbands], sigimagelist[i]);
                strcpy(cmasklist[*nbands], cmasklist[i]);
                strcpy(svblist[*nbands], svblist[i]);
            }
            bandlist[*nbands] = i;
            (*nbands)++;
        }
    }

    if (*nbands == 0)
    {
        printf("MDET: No input images specified.\n");
        return 2;
    }
    if (*backwindow == 0)
    {
        printf("MDET: I need a BACKWINDOW size [coadd pixels] (float).\n");
        return 2;
    }
    if (*threshold == 0)
    {
        printf("MDET: I need a THRESHOLD [sigmas] (float).\n");
        return 2;
    }
    if (koutlist == -1)
    {
        printf("MDET: I need an OUTLIST (character string).\n");
        return 2;
    }

    return 0;
}