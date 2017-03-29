#pragma once

#include "__base.h"

int weedout_mf(float* x, float* y, float* sc, int nsources, float rcrit, float pixel, BOOL* notseen,
    int* satmask, float* rsatset, int nx, int ny, int nbands, int isingleband);