#pragma once

#define MAX_BANDS 4

#define BLANK_PIXEL -999.0
#define MAX_SAT 100

#define MAX_SOURCES 1000000
#define MAX_ANN 10000
/*
struct Band
{
    char imagepath[MAX_PATH], sigimagepath[MAX_PATH], cmaskpath[MAX_PATH];
    char svbpath[MAX_PATH];
    float fwhm, focalpix, satradius;
};
*/

struct Image
{
    int nx, ny;
    float ra, dec, cdelt1, cdelt2, theta;
    float* pixels;
};

struct Mask
{
    int nx, ny;
    float ra, dec, cdelt1, cdelt2, theta;
    int* pixels;
};

struct Source
{
    float x, y, value;
};

int mdet(int argc, char** argv);