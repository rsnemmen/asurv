ASURV: Astronomy SURVival analysis (Rev. 1.3)
==============================================

This shell archive contains source code and documentation for ASURV,
Astronomy SURVival analysis (Rev. 1.3).  ASURV implements a suite of
statistical methods for the analysis of censored data; i.e. data
which are known to lie above or below some limit.   It was written
specifically to treat left-censoring arising in observational astronomy
when objects are observed but sometimes not detected due to sensitivity 
limits.  However, the methods can be useful to researchers in other 
disciplines, as the code includes techniques that are often omitted 
from commercial survival analysis packages. 

ASURV computes: 

- the maximum-likelihood Kaplan-Meier estimator
- several univariate two-sample tests (Gehan, Peto-Peto, Peto-Prentice)
- three bivariate correlation coefficients (Cox regression, generalized Kendall's tau and Spearman's rho)
- three linear regressions (EM algorithm assuming normal residuals, Buckley-James line, Schmitt line).  

The user 
is strongly encouraged to read the manual and references below if they
are unfamiliar with these methods.  The program is stand-alone and does
not call any specialized library.

# Files included

This archive contains: 

- this `README` file
- `asurv.f` with 59 surbroutines in FORTRAN 77
- `asurv.doc` and `asurv.tex` with the Users' Manual in ASCII and LaTeX respectively
- `asurv.etc` with test files, subroutine list and Bug Report form.  

These files were retrieved from the [Center for AstroStatistics at PennState website](https://astrostatistics.psu.edu/statcodes/asurv).  

# Compilation

The source code can be compiled with
a statement like 

    gfortran -o asurv asurv.f  

~~Operation has been verified for a variety of computer platforms including UNIX, VMS, VM and DOS.~~ Please verify if this works in your favorite platform and submit a PR if not with improvements.

# Authors

ASURV was written between 1987 and 1992 by Drs. Takashi Isobe (Center
for Space Research, MIT), Michael LaValley (formerly at Dept. of 
Statistics, Penn State), and Eric Feigelson.  Code development was 
supported by several NASA grants.  Questions and problems should be 
addressed to:  Eric Feigelson, Dept. of Astronomy & Astrophysics, 
Pennsylvania State University, University Park PA USA, FAX 814-863-3399, 
Email edf@astro.psu.edu, WWW http://www.astro.psu.edu/users/edf).   

IMPORTANT: The authors grant researchers and students permission to
use and copy ASURV code and associated material for non-commercial
purposes.  We request that publications resulting from its use cite
one of the references below.  This software is provided `as is' without
any expressed or implied warranty.  

References:

- Feigelson, E. D. and Nelson, P. I. Statistical Methods for Astronomical Data with Upper Limits: I. Univariate Distributions, Astrophyscal Journal 293, 192-206, 1985
- Isobe, T., Feigelson, E. D., and Nelson, P. I. Statistical Methods for Astronomical Data with Upper Limits: II. Correlation and Regression, Astrophysical Journal, 306, 490-507, 1986
- LaValley, M., Isobe, T. and Feigelson, E.D. ``ASURV'', Bulletin American Astronomical Society (Software Reports),  22, 917-918, 1990

Revisions:
Rev. 0 (1987-1990)  Incomplete and obsolete.
Rev. 1 (1992-present) Essentially identical versions with minor bugs corrected. 
 
README for asurv.shar written by Eric feigelson (Sept. 1996)

NOTE: Some users have encountered difficulty compiling ASURV on
Linux and MacOS systems.  This problem may be solved by renaming subroutines STAT and UNPACK throughout the code.  



