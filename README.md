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

Use `make` to build the standalone executable:

    make

The Makefile defaults to `gfortran` and adds the legacy-compatibility
flags needed by modern compilers. The equivalent manual command is:

    gfortran -std=legacy -fallow-argument-mismatch -o asurv asurv.f

Useful convenience targets:

- `make help` shows the available targets
- `make clean` removes the built executable
- `make distclean` also removes common compiler scratch files
- `make test` runs the bundled smoke and exact-output regression tests

# Validation

Bundled smoke tests based on the examples in `asurv.etc` are available as:

    make smoke
    make smoke-gal1
    make smoke-gal2
    make smoke-gal3

Exact-output regression tests using the embedded sample outputs in `asurv.etc`
are also available:

    make regression
    make regression-gal1
    make regression-gal2
    make regression-gal3

# Usage

[![asciicast](https://asciinema.org/a/265869.svg)](https://asciinema.org/a/265869)

An experimental Python wrapper is also available as `asurv_cli.py`. It keeps
the legacy Fortran executable as the backend, but replaces command-file authoring
and menu navigation with subcommands. It currently expects the existing
legacy-formatted ASURV data files as input.

Examples:

    python3 asurv_cli.py km \
      --data-file gal1.dat \
      --title "IR Luminosities of Galaxies" \
      --nvar 1 \
      --column 1 \
      --variable-name IR \
      --full-km \
      --diff-start 25 \
      --diff-bins 5 \
      --diff-size 2 \
      --print-data

    python3 asurv_cli.py twost \
      --data-file gal2.dat \
      --title "IR Luminosities of Galaxies" \
      --nvar 1 \
      --column 1 \
      --variable-name IR \
      --group-id 0 --group-id 1 --group-id 2 \
      --first-group 0 \
      --second-group 1 \
      --first-group-name Normal \
      --second-group-name Starburst \
      --full-km

    python3 asurv_cli.py bivar \
      --data-file gal3.dat \
      --title "Optical-Infrared Relation" \
      --n-independent 1 \
      --x-name Optical \
      --y-name Infrared \
      --method cox \
      --method em \
      --em-tolerance 1.0e-5 \
      --em-initial-coefficients 0 0 0 \
      --em-max-iterations 50

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
