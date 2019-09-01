NOTE: Some users have encountered difficulty compiling ASURV on
Linux and MacOS systems.  This problem may be solved by renaming
subroutines STAT and UNPACK throughout the code.  

XREADME for ASURV
X
XThis shell archive contains source code and documentation for ASURV,
XAstronomy SURVival analysis (Rev. 1.3).  ASURV implements a suite of
Xstatistical methods for the analysis of censored data; i.e. data
Xwhich are known to lie above or below some limit.   It was written
Xspecifically to treat left-censoring arising in observational astronomy
Xwhen objects are observed but sometimes not detected due to sensitivity 
Xlimits.  However, the methods can be useful to researchers in other 
Xdisciplines, as the code includes techniques that are often omitted 
Xfrom commercial survival analysis packages. 
X
XASURV computes: the maximum-likelihood Kaplan-Meier estimator; several 
Xunivariate two-sample tests (Gehan, Peto-Peto, Peto-Prentice); three 
Xbivariate correlation coefficients (Cox regression, generalized Kendall's
Xtau and Spearman's rho); and three linear regressions (EM algorithm
Xassuming normal residuals, Buckley-James line, Schmitt line).  The user 
Xis strongly encouraged to read the manual and references below if they
Xare unfamiliar with these methods.  The program is stand-alone and does
Xnot call any specialized library.
X
XThis archive contains: this README file; asurv_code.f with 59
Xsurbroutines in FORTRAN 77; asurv.doc and asurv.tex with the Users'
XManual in ASCII and LaTeX respectively; and asurv.etc with test files,
Xsubroutine list and Bug Report form.  It is distributed via the World
XWide Web (ftp://www.astro.psu.edu/users/edf/asurv.shar) or by email
Xrequest to code@stat.psu.edu.  The archive is unpacked using a UNIX
Xcommand like `sh asurv.shar'.  The source code can be compiled with
Xa statement like `f77 -f -o asurv asurv_code.f'.  Operation has been
Xverified for a variety of computer platforms including UNIX, VMS, VM 
Xand DOS.
X
XASURV was written between 1987 and 1992 by Drs. Takashi Isobe (Center
Xfor Space Research, MIT), Michael LaValley (formerly at Dept. of 
XStatistics, Penn State), and Eric Feigelson.  Code development was 
Xsupported by several NASA grants.  Questions and problems should be 
Xaddressed to:  Eric Feigelson, Dept. of Astronomy & Astrophysics, 
XPennsylvania State University, University Park PA USA, FAX 814-863-3399, 
XEmail edf@astro.psu.edu, WWW http://www.astro.psu.edu/users/edf).   
X
XIMPORTANT: The authors grant researchers and students permission to
Xuse and copy ASURV code and associated material for non-commercial
Xpurposes.  We request that publications resulting from its use cite
Xone of the references below.  This software is provided `as is' without
Xany expressed or implied warranty.  
X
XReferences:
XFeigelson, E. D. and Nelson, P. I. ``Statistical Methods for
X    Astronomical Data with Upper Limits: I. Univariate Distributions',
X    Astrophyscal Journal 293, 192-206, 1985.
XIsobe, T., Feigelson, E. D., and Nelson, P. I. ``Statistical Methods
X    for Astronomical Data with Upper Limits: II. Correlation and Regression',
X    Astrophysical Journal, 306, 490-507, 1986.
XLaValley, M., Isobe, T. and Feigelson, E.D. ``ASURV'', Bulletin
X    Amercan Astronomical Society (Software Reports),  22, 917-918, 1990.
X
XRevisions:
XRev. 0 (1987-1990)  Incomplete and obsolete.
XRev. 1 (1992-present) Essentially identical versions with minor bugs
X    corrected. 
X 
XREADME for asurv.shar written by Eric feigelson (Sept. 1996)
X