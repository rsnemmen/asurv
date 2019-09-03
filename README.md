ASURV: Astronomy SURVival analysis 
=======================================

This contains source code and documentation for ASURV,
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
- `manual.txt` and `manual.tex` with the Users' Manual in ASCII and LaTeX respectively
- `asurv.etc` with test files, subroutine list and Bug Report form.  

These files were retrieved from the [Center for AstroStatistics at PennState website](https://astrostatistics.psu.edu/statcodes/asurv) and updated to fit modern architectures.

# Compilation

The source code can be compiled with
a statement like 

    gfortran -o asurv asurv.f  

# Usage

[![asciicast](https://asciinema.org/a/265869.svg)](https://asciinema.org/a/265869)

# References

If you use ASURV, you are morally obliged to cite the following papers written by the authors of the code:

- Feigelson, E. D. and Nelson, P. I. Statistical Methods for Astronomical Data with Upper Limits: I. Univariate Distributions, Astrophyscal Journal 293, 192-206, 1985
- Isobe, T., Feigelson, E. D., and Nelson, P. I. Statistical Methods for Astronomical Data with Upper Limits: II. Correlation and Regression, Astrophysical Journal, 306, 490-507, 1986
- LaValley, M., Isobe, T. and Feigelson, E.D. ASURV, Bulletin American Astronomical Society (Software Reports),  22, 917-918, 1990

# Revisions

- Rev. 0 (1987-1990)  Incomplete and obsolete.
- Rev. 1 (1992-present) Essentially identical versions with minor bugs corrected. 
- Rev. 3 (9/2019) Code made available on github and compilation bugs fixed
 
 README written by Eric feigelson (Sep. 1996) and revised by [Rodrigo Nemmen](https://rodrigonemmen.com) (Sep. 2019).
