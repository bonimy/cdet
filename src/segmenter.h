#pragma once

void segmenter(int* msk, int nx, int ny, int iblo, int ibhi, int jblo, int jbhi, int isearchrad,
    int *ixsat, int *iysat, float* satradius, int* nsatpix);
int fillneighbors(int* msk, int* objmsk, int nx, int ny, int iblo, int ibhi, int jblo, int jbhi,
	int sz, int x, int y, int oindex);