#pragma once

#include "__base.h"

void wimage(int nx, int ny, float* larray, char* fin, char* fout);
void deletefile(char* filename, int* status, BOOL zexist, BOOL erase);