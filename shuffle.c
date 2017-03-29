//Rearrange 1-d vector LARRAY as a 2-d ARRAY.

#include "__base.h"
#include "shuffle.h"

void shuffle(float* larray, int nx, int ny, float* array)
{
    memcpy(array, larray, nx * ny * sizeof(float));
}