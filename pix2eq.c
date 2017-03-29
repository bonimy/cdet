//Transformation from pixel coordinates to equatorial(RA, Dec) based on sine
//projection.
//
//Input parameters :
//     i, j = pixel location
//     ra0, dec0 = reference position[deg]
//     iref, jref = reference pixel
//     pixel = pixel size[deg]
//     theta = rotation angle(CROTA2)
//
//Output parameters :
//     ra, dec = equatorial position[deg]

#include "__base.h"
#include "pix2eq.h"

void pix2eq(double i, double j, double ra0, double dec0, double iref, double jref,
    double pixel, double theta, double* ra, double* dec)
{
    const double dtor = 0.01745329251994;

    i -= iref;
    j -= jref;

    ra0 *= dtor;
    dec0 *= dtor;
    theta *= dtor;
    i *= dtor * pixel;
    j *= dtor * pixel;

    double x = i * cos(theta) + j * sin(theta);
    double y = i * sin(theta) - j * cos(theta);

    double D = asin(sqrt(x * x + y * y));
    double B = atan2(-x, y);
    
    x = sin(dec0) * sin(D) * cos(B) + cos(dec0)*cos(D);
    y = sin(D) * sin(B);

    ra0 += atan2(y, x);
    dec0 = asin(sin(dec0) * cos(D) - cos(dec0) * sin(D) * cos(B));

    *ra = ra0 / dtor;
    *dec = dec0 / dtor;
}