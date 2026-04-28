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

# Python Interface

## Install

Build the Fortran backend first, then install the Python package in the same directory:

```bash
make
pip install -e ".[plot,dev]"
```

## Importing

```python
import pandas as pd
import asurv

# Kaplan-Meier
df = pd.read_fwf("gal1.dat", widths=[4, 10], header=None, names=["ind", "lum"])
df["upper_limit"] = df["ind"] < 0

result = asurv.km(df, value="lum", upper_limit="upper_limit",
                  title="IR Luminosities of Galaxies",
                  differential=(25.0, 5, 2.0))

print(result.summary())
print(result.table)        # pandas DataFrame: from, to, S, error
print(result.mean)         # 27.767
result.plot()              # survival curve with censoring ticks

# Two-sample tests
df2 = pd.read_fwf("gal2.dat", widths=[4, 4, 10], header=None,
                   names=["group", "ind", "lum"])
df2["upper_limit"] = df2["ind"] < 0

r2 = asurv.twosample(df2, value="lum", group="group", upper_limit="upper_limit",
                     groups=(0, 1), title="Normal vs Starburst")

print(r2.summary())                      # table of all 5 test statistics
print(r2.tests["logrank"].probability)   # p-value
r2.plot(labels=("Normal", "Starburst"))  # overlaid KM curves

# Bivariate correlation & regression
df3 = pd.read_fwf("gal3.dat", widths=[4, 10, 10], header=None,
                   names=["ind", "optical", "infrared"])
df3["upper_limit"] = df3["ind"] < 0

r3 = asurv.bivar(df3, x="optical", y="infrared",
                 methods=["cox", "em"], y_upper="upper_limit",
                 title="Optical-Infrared Relation",
                 x_name="Optical", y_name="Infrared")

print(r3.summary())              # Cox chi², EM intercept/slope/sigma
print(r3.cox.chi2)               # 18.458
print(r3.em.slopes)              # [(1.0519, 0.0793)]
r3.plot(df3, x="optical", y="infrared", y_upper="upper_limit")
```

## CLI

The `asurv-py` command is installed automatically:

```bash
asurv-py km --data-file gal1.dat --title "IR Luminosities" --full-km \
            --diff-start 25 --diff-bins 5 --diff-size 2

asurv-py twost --data-file gal2.dat --title "Two-Sample Test" \
               --first-group-name Normal --second-group-name Starburst

asurv-py bivar --data-file gal3.dat --title "Optical-IR" \
               --x-name Optical --y-name Infrared \
               --method cox --method em --em-initial-coefficients 0 0 0
```

The legacy `python3 asurv_cli.py ...` invocations continue to work via a backward-compatibility shim.

## Running tests

```bash
make test          # Fortran smoke + regression tests (unchanged)
pytest             # Python unit and integration tests
pytest -m "not integration"  # fast offline tests only (no binary needed)
```

## Batch mode (direct)

The Fortran binary can still be called directly in batch mode:

    ./asurv KM km.com
    ./asurv TWOST ts.com
    ./asurv BIVAR bv.com

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
