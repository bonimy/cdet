# CDET - Multi-wavelength source detection software ([MDET](https://github.com/catwise/mdet) written in C)
For information on what the MDET algorithm is, take look at [Optimal Multiwavelength Source Detection: Experience Gained
from the WISE Mission (K. Marsh and T. Jarrett)](MarshJarrett2012.pdf).

# History
The MDET module was first written in 2008 for a part of the WISE catalog software pipeline at IPAC. Because the WISE computer cluster is occupied with ongoing satellite operations, we need to get the software running on independent systems not on the WISE cluster.
The last iteration of MDET was written in Fortran. CDET is a C-port written by Nelson Garcia.

# Building and Compiling
CDET is written in the C programming language. CDET's only dependency (besides the standard libraries) is the [FITSIO Subroutine Library](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html). The CFITSIO (C Libraries) are written to always be backwards compatible, so simply download the latest version and install it normally.

The build commands we use for CDET are

`gcc *.c -o cdet -m64 -O3 -std=c99 -L../cfitsio -lcfitsio`

Replace `-L../cfitsio` with the directory where you are keeping the CFITSIO library. If this library is already registered in your libpath, then you do not need it, just `-lcfitsio`. Be sure to run this command in the same directory as the CDET source files. `*.c` tell GCC (The GNU C compiler) to compile all source code file for CDET.


# Using MDET
[ToDo: Document usage]

# License
MDET is written under the [GNU General Public License v3.0](LICENSE).

# Code of Conduct
Contributors to MDET are required to follow the [Contributor Covenant Code of Conduct](code_of_conduct.md)
