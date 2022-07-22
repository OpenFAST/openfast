                                 FITPACK

FITPACK is a collection of FORTRAN programs for CURVE and SURFACE FITTING with
SPLINES and TENSOR PRODUCT SPLINES. 

Features included are automatic knot selection, error smoothing and data 
reduction. Besides the set of data points, the user merely has to provide a
single parameter, called the smoothing factor, to control the tradeoff between
closeness of fit and smoothness of fit. Confidence limits for this parameter
are available if the statistical errors on the data can be estimated.

In addition to a standard routine for curve fitting over an interval and a  
routine for fitting a bivariate spline function to a set of scattered data on a
rectangular domain, there are special programs for smoothing periodic functions,
parametric curves and surfaces, functions with convexity constraints. data over
the sphere and other non-rectangular approximation domains. Very efficient
surface fitting routines are provided for applications with data given either
on a rectangular, polar or spherical grid. Routines are included for evaluation,
differentiation and integration of the spline approximations. FITPACK also
contains an insertion algorithm and routines for calculating zeros and Fourier 
coefficients of cubic splines.

The subdirectory ex contains a test program for each of the 29
main routines, e.g. the file mncurf.f contains a sample program for curfit,
mnperc.f a test program for percur.f, etc. Some of the test programs need some
input data. These can be found in the files with initials da, e.g. the program
in mnpasu.f needs the data which are stored in dapasu, etc.

The user-level routines should be fairly well documented. Besides a parameter
description, there is a pointer to the lower level FITPACK routines, needed by
this particular routine. This may be useful in case only a restricted number
of fitting routines would be needed. A user guide of approximately 300 pages
can be ordered by sending a check of 25 US dollars to
		Prof. P. Dierckx
		Department of Computer Science
		K.U.Leuven
		Celestijnenlaan 200 A
		B 3001, Heverlee, Belgium

The mathematical fundamentals of FITPACK are described in the book:
   Curve and Surface Fitting with Splines, by P. Dierckx,
   Monographs on Numerical Analysis, Oxford University Press, 1993,
   (ISBN 0-19-853441-8, 286 pages, 56 line figures)
In addition, potential FITPACK users will find there a lot of practical advice,
illustrated with many examples, both purely academic as well as taken from
real live. Comments on the software or the book may be addressed to the author
via email at the following address: Paul.Dierckx@cs.kuleuven.ac.be
                               
Finally, it should be noted that this software has no connection with
the FITPACK software package of Alan Cline.
