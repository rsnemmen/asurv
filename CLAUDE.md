# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make          # build the asurv executable
make clean    # remove executable
make distclean  # remove executable and compiler scratch files (*.o, *.mod, fort.*)
```

The Makefile uses `gfortran -std=legacy -fallow-argument-mismatch`. Equivalent manual command:

```bash
gfortran -std=legacy -fallow-argument-mismatch -o asurv asurv.f
```

## Tests

```bash
make test              # run all smoke + regression tests
make smoke             # smoke tests (output sanity check only)
make regression        # exact-output diff against expected in asurv.etc
make smoke-gal1        # Kaplan-Meier sample
make smoke-gal2        # two-sample test sample
make smoke-gal3        # bivariate sample
make regression-gal1   # Kaplan-Meier exact-output
make regression-gal2   # two-sample exact-output
make regression-gal3   # bivariate exact-output
```

The Makefile extracts data, command, and expected-output sections from `asurv.etc` into a temp directory, runs the executable via piped stdin menu choices, then diffs the output.

## Architecture

The entire Fortran application lives in `asurv.f` — one fixed-form Fortran 77 source file (~59 subroutines). There is no multi-file build; the old subroutine list in `asurv.etc` is historical only.

**Top-level dispatch:**
- `MAIN` shows the top-level menu and dispatches to `UNIVAR` (menu choice 1) or `BIVAR` (menu choice 2).
- `UNIVAR` handles Kaplan-Meier estimation and two-sample tests. It reads data via `DATAIN`, calls `KMESTM` for Kaplan-Meier and `TWOST` for two-sample statistics.
- `BIVAR` handles correlation and regression. It reads data via `DATREG`, then dispatches to `COXREG`, `BHK`, `SPRMAN`, `EM`, `BJ`, or `TWOKM`. Method-specific prompts are in helper routines `R3`–`R6`.

**Batch modes** (bypass the interactive menu):

```bash
./asurv KM km.com
./asurv TWOST ts.com
./asurv BIVAR bv.com
```

**Python CLI wrapper** (`asurv_cli.py`): an experimental Python CLI that generates command files and drives the Fortran executable in batch mode. Subcommands: `km`, `twost`, `bivar`.

**Automating runs:** prefer command files plus piped top-level menu choices over scripting every interactive prompt. The authoritative command-file formats are defined in `UNIVAR`, `BIVAR`, and `R3`–`R6` in `asurv.f`; worked examples are in `manual.txt` and `asurv.etc`.

## Key Fortran Conventions

- Fixed-form Fortran 77: column-1 `C` comments, `+` continuation lines, `IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)` throughout numerical routines.
- Hard size limits declared as `PARAMETER(MVAR=4, NDAT=500, IBIN=50)` in **both** `UNIVAR` and `BIVAR`. If raising limits, update both.
- Input is fixed-width: `DATAIN` uses `FORMAT(10(I4,F10.3))` for KM data and `FORMAT(I4,10(I4,F10.3))` for two-sample data; `DATREG` uses `FORMAT(I4,11F10.3)` for bivariate data. If `MVAR` grows past 10, widen these formats per Appendix A7 of the manual. If data needs more than three decimal places, update the relevant formats and reduce `CONST` in `CENS` so it retains at least two more decimal places than the data.
- `CHARACTER*9` is a hard limit for filenames, variable names, group labels, command files, and output files.
- Do not remove the upper-limit sign normalization in `AARRAY` and `XVAR` — these routines flip upper-limit data via `ISIGN` so downstream statistics can reuse lower-limit logic.
- Work arrays are allocated in `UNIVAR`/`BIVAR` and passed through long argument lists rather than local allocation.
- `TWOST` always computes all tests (Gehan permutation, Gehan hypergeometric, logrank, Peto-Peto, Peto-Prentice) — this is intentional, not redundant.

## Sample Data

`asurv.etc` contains embedded sections delimited by `****` headers:
- `gal1.dat` / `gal1.com` / `gal1.out` — Kaplan-Meier example (IR Luminosities)
- `gal2.dat` / `gal2.com` / `gal2.out` — two-sample test example
- `gal3.dat` / `gal3.com` / `gal3.out` — bivariate correlation/regression example
