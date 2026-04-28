# ASURV: Astronomy SURVival analysis

ASURV is a survival-analysis package for censored data, especially upper limits
common in observational astronomy. This repository preserves the original
Fortran implementation and now ships a maintained Python interface that wraps
the backend for programmatic use and a modern CLI.

**Most users should start with the Python package.** The `asurv` executable
remains the numerical backend and the authoritative implementation of the
statistics, while the Python layer handles data preparation, backend execution,
and parsing the text output into structured results.

## What ASURV computes

ASURV provides:

- the maximum-likelihood Kaplan-Meier estimator
- several univariate two-sample tests (Gehan, logrank, Peto-Peto, Peto-Prentice)
- bivariate correlation methods (Cox regression, generalized Kendall's tau,
  Spearman's rho)
- linear regressions for censored data (EM algorithm, Buckley-James, Schmitt)

If you are not already familiar with these methods, read `manual.txt` and the
papers cited below before using the results scientifically.

## Interfaces

| Interface | Best for | Entry points |
| --- | --- | --- |
| Python API | notebooks, scripts, reproducible analyses | `asurv.km()`, `asurv.twosample()`, `asurv.bivar()` |
| CLI | shell workflows without writing Python code | `asurv-py km`, `asurv-py twost`, `asurv-py bivar` |
| Direct backend | legacy workflows, command files, regression testing | `./asurv KM ...`, `./asurv TWOST ...`, `./asurv BIVAR ...` |

The Python package wraps the same Fortran backend used by the direct interface.
For numerical correctness, treat the Fortran executable as the source of truth.

## Installation

### Prerequisites

- Python 3.9+
- `gfortran` or another compiler compatible with the Makefile

### Build the backend

```bash
make
```

The Makefile uses the legacy-compatibility flags required by modern `gfortran`
versions. The equivalent manual command is:

```bash
gfortran -std=legacy -fallow-argument-mismatch -o asurv asurv.f
```

### Install the Python package

For normal use:

```bash
pip install -e .
```

If you want plotting support:

```bash
pip install -e ".[plot]"
```

If you are contributing or running the full developer workflow:

```bash
pip install -e ".[plot,dev]"
```

The Python package requires the Fortran binary at runtime. By default it looks
for:

1. `ASURV_BIN`
2. `./asurv` in the current working directory
3. `asurv` on `PATH`

## Python quick start

The package exposes three high-level entry points:

- `asurv.km()` for Kaplan-Meier estimation
- `asurv.twosample()` for univariate two-sample tests
- `asurv.bivar()` for bivariate correlation and regression

All three return structured result objects with `.summary()`, `.raw_output`,
and method-specific parsed fields. Plotting helpers are available when
installed with `.[plot]`.

### Kaplan-Meier

```python
import pandas as pd
import asurv

df = pd.read_fwf("gal1.dat", widths=[4, 10], header=None, names=["ind", "lum"])
df["upper_limit"] = df["ind"] < 0

result = asurv.km(
    df,
    value="lum",
    upper_limit="upper_limit",
    title="IR Luminosities of Galaxies",
    differential=(25.0, 5, 2.0),
)

print(result.summary())
print(result.table)
print(result.mean)
result.plot()
```

### Two-sample tests

```python
import pandas as pd
import asurv

df = pd.read_fwf(
    "gal2.dat",
    widths=[4, 4, 10],
    header=None,
    names=["group", "ind", "lum"],
)
df["upper_limit"] = df["ind"] < 0

result = asurv.twosample(
    df,
    value="lum",
    group="group",
    upper_limit="upper_limit",
    groups=(0, 1),
    title="Normal vs Starburst",
)

print(result.summary())
print(result.tests["logrank"].probability)
result.plot(labels=("Normal", "Starburst"))
```

### Bivariate correlation and regression

```python
import pandas as pd
import asurv

df = pd.read_fwf(
    "gal3.dat",
    widths=[4, 10, 10],
    header=None,
    names=["ind", "optical", "infrared"],
)
df["upper_limit"] = df["ind"] < 0

result = asurv.bivar(
    df,
    x="optical",
    y="infrared",
    methods=["cox", "em"],
    y_upper="upper_limit",
    title="Optical-Infrared Relation",
    x_name="Optical",
    y_name="Infrared",
)

print(result.summary())
print(result.cox.chi2)
print(result.em.slopes)
result.plot(df, x="optical", y="infrared", y_upper="upper_limit")
```

## CLI

Installing the package exposes `asurv-py`, which mirrors the Python API:

```bash
asurv-py km --data-file gal1.dat --title "IR Luminosities" --full-km \
            --diff-start 25 --diff-bins 5 --diff-size 2

asurv-py twost --data-file gal2.dat --title "Two-Sample Test" \
               --first-group-name Normal --second-group-name Starburst

asurv-py bivar --data-file gal3.dat --title "Optical-IR" \
               --x-name Optical --y-name Infrared \
               --method cox --method em --em-initial-coefficients 0 0 0
```

Legacy invocations through `python3 asurv_cli.py ...` still work via a
backward-compatibility shim, but new usage should prefer `asurv-py`.

## Validation and tests

### Build and maintenance targets

- `make help` shows available Make targets
- `make clean` removes the built executable
- `make distclean` also removes common compiler scratch files

### Fortran regression suite

The bundled sample cases in `asurv.etc` drive the backend validation suite:

```bash
make test
make smoke
make regression
make smoke-gal1
make smoke-gal2
make smoke-gal3
make regression-gal1
make regression-gal2
make regression-gal3
```

The Fortran suite is the authoritative check for numerical correctness.

### Python test suite

```bash
pytest
pytest -m "not integration"
```

The full pytest suite exercises the Python wrapper, parser, and integration with
the Fortran backend. The offline subset skips tests that require a built
`./asurv` binary.

## Direct Fortran usage

The standalone executable still supports direct batch execution:

```bash
./asurv KM km.com
./asurv TWOST ts.com
./asurv BIVAR bv.com
```

This corresponds to the traditional ASURV workflows:

- `UNIVAR` for Kaplan-Meier estimation and two-sample tests
- `BIVAR` for correlation and regression

If you are working directly with command files, fixed-width input files, or
legacy output reports, use `manual.txt` and `manual.tex` as the primary
reference.

## Repository layout

| Path | Purpose |
| --- | --- |
| `asurv.f` | original fixed-form Fortran implementation |
| `asurv.etc` | bundled sample data, command files, expected outputs |
| `manual.txt`, `manual.tex` | user manual |
| `src/asurv/` | Python wrapper, CLI, backend launcher, parsers, plotting helpers |
| `tests/` | pytest suite for wrapper, parser, validation, and integration coverage |
| `asurv_cli.py` | backward-compatibility shim to `src/asurv/cli.py` |
| `pyproject.toml` | Python package metadata and console-script definition |

These files were originally retrieved from the Center for Astrostatistics at
Penn State and later updated to build on modern systems and support the Python
wrapper.

## References

If you use ASURV, please cite:

- Feigelson, E. D. and Nelson, P. I. Statistical Methods for Astronomical Data
  with Upper Limits: I. Univariate Distributions, Astrophysical Journal 293,
  192-206, 1985
- Isobe, T., Feigelson, E. D., and Nelson, P. I. Statistical Methods for
  Astronomical Data with Upper Limits: II. Correlation and Regression,
  Astrophysical Journal 306, 490-507, 1986
- LaValley, M., Isobe, T. and Feigelson, E. D. ASURV, Bulletin of the American
  Astronomical Society (Software Reports), 22, 917-918, 1990

## Revisions

- Rev. 0 (1987-1990): incomplete and obsolete
- Rev. 1 (1992-present): essentially identical versions with minor bugs corrected
- Rev. 3 (2019): code made available on GitHub and compilation bugs fixed

README originally written by Eric Feigelson (Sep. 1996) and later revised by
[Rodrigo Nemmen](https://rodrigonemmen.com).
