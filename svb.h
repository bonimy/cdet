#pragma once

void svb(float* a, float* unca, int nx, int ny, int width, int iskip, float* b, float* csig);

int fillbox(float* a, float* unca, int nx, int ny, float* box, int width, int x, int y, int iskip, float* sum);
void interpolate(float* b, float* csig, float* bwt, int nx, int ny, float* gauss, int ik, float vmed, float vcon, int x, int y);