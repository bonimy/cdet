#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI 3.14159265358979323846
#endif
#endif // !PI

#define RAD_PER_DEG (M_PI/180)
#define DEG_PER_RAD (180/M_PI)

#ifndef MAX_PATH
#ifdef _MAX_PATH
#define MAX_PATH 260
#else // !_MAX_PATH
#define MAX_PATH 260
#endif // _MAX_PATH
#endif // !MAX_PATH

#ifndef BOOL
#define BOOL int
#endif // !BOOL

#ifndef TRUE
#define TRUE -1
#endif // !TRUE

#ifndef FALSE
#define FALSE 0
#endif // !FALSE

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif // !max

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif // !min