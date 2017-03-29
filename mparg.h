#pragma once

#include "__base.h"
#include "mdet.h"

int mparg(char** argv, int argc, int* nbands, int* backwindow, float* threshold,
    float fwhm[MAX_BANDS], float focalpix[MAX_BANDS], BOOL bandlist[MAX_BANDS],
    char imagelist[MAX_BANDS][MAX_PATH], char cmasklist[MAX_BANDS][MAX_PATH],
    char svblist[MAX_BANDS][MAX_PATH], char sigimagelist[MAX_BANDS][MAX_PATH],
    char outlist[MAX_PATH], char matchout[MAX_PATH], BOOL* confused, BOOL* multiframe, int* sigfigs);