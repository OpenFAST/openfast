

ROOTS OF QUARTIC/CUBIC/QUADRATIC POLYNOMIALS                 /quartic/readme.txt
   WITH REAL COEFFICIENTS
                              

The files for this program are in the directory \quartic on the CD-ROM and in the
archive file quartic.zip that may be downloaded from the PDAS web site.
  readme.txt       this file of general description
  quartic.f90      the complete source codetosh (Intel)


The reference documents for this program may be accessed
from the web page http://www.pdas.com/quarticrefs.html. 

2020-03-05 bjj: see http://www.pdas.com/quartic.html and http://www.pdas.com/quarticdownload.html

To compile this program for your computer, use the command
   gfortran  quartic.f90 -o quartic.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.


Many problems in science and engineering lead to polynomial equations
and the desired physical quantities must be found by solving for the
zeroes of the equation. Most books on Numerical Computing or Engineering
Mathematics show examples of code for making these calculations.
One must be careful with roundoff and overflow when making these
calculations and the textbook examples frequently do not incorporate
these "robustness" features.

There is a collection of software from the U.S. Naval Surface Weapons
Center that has been widely distributed and checked for accuracy.
Among this collection is a very nice coding of the solution of zeroes 
of polynomial equations with real coefficients up to quartic order. 
I noted that Alan Miller of CSIRO had updated the code to comply with 
modern Fortran, using the Essential Lahey Fortran 90 compiler, which 
enforces very strict standards of program structure and syntax. 
I have added a simple front-end that allows you to solve quartics as 
a stand-alone program. Of course, you may want to extract the subroutine 
for inclusion in your own code.

This program computes the solution to the polynomial equation
  a*x**4 + b*x**3 + c*x**2 + d*x + e = 0
with real coefficients.

The program asks for the coefficients of each term of the polynomial
to be solved. If you are solving a cubic, answer zero to the question 
"What is the coefficient of the quartic term?"

If you want to take the program apart to use in your own programs, you
can pick up each of the individual subroutines or go for the general
subroutine SolvePolynomial. It should be pretty obvious. If not, send
me mail and I will try to clear it up.

The routine quartic.f90 has been checked and will compile properly with 
either of the free Fortran compilers, gfortran or g95.

There is a treasure house of wonderful Fortran code at Alan Miller's web site
   at http://users.bigpond.net.au/amiller/
There seems to be a problem with this web site, but there is a mirror site at
   http://jblevins.org/mirror/amiller/
that should work for now.

