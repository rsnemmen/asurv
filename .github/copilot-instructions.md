# ASURV Copilot Instructions

## Build and verification

- Build the standalone executable from the consolidated Fortran 77 source:
  ```bash
  gfortran -o asurv asurv.f
  ```
- There is no dedicated lint target or automated test suite in this repository. Use the bundled examples in `asurv.etc` as smoke tests.
- Single smoke test for the Kaplan-Meier path (`gal1.dat` / `gal1.com` embedded in `asurv.etc`):
  ```bash
  repo_root=$(pwd)
  tmpdir=$(mktemp -d)
  trap 'rm -rf "$tmpdir"' EXIT
  awk 'NR>=66&&NR<=71{print}' asurv.etc > "$tmpdir/gal1.dat"
  awk 'NR>=75&&NR<=86{print}' asurv.etc > "$tmpdir/gal1.com"
  (
    cd "$tmpdir" &&
    printf '1\n1\ny\ngal1.com\nn\n' | "$repo_root"/asurv >/dev/null &&
    grep -q 'KAPLAN-MEIER ESTIMATOR' gal1.out
  )
  ```
- `asurv.etc` also contains reusable TWOST (`gal2.*`) and BIVAR (`gal3.*`) sample data, command files, and expected output.

## High-level architecture

- `asurv.f` is the whole application: one fixed-form Fortran source file containing the main menu plus all statistical subroutines. The old multi-file build list is preserved only as historical reference inside `asurv.etc`; the maintained build path in this repo is the single-file `gfortran -o asurv asurv.f` flow from `README.md`.
- The main program shows the top-level menu, then dispatches to `UNIVAR` for univariate workflows or `BIVAR` for correlation/regression workflows.
- `UNIVAR` owns both Kaplan-Meier estimation and two-sample testing. It reads fixed-width data through `DATAIN`, gathers menu/command-file input itself, calls `KMESTM` for Kaplan-Meier output, and calls `TWOST` for two-sample statistics before optionally emitting per-group Kaplan-Meier summaries.
- `BIVAR` reads bivariate data through `DATREG`, decides which methods are legal from the censoring pattern and number of independent variables, then dispatches to `COXREG`, `BHK`, `SPRMAN`, `EM`, `BJ`, or `TWOKM`. Method-specific extra prompts are isolated in `R3`, `R4`, `R5`, and `R6`.
- Most algorithm routines work on fixed-size work arrays allocated in `UNIVAR` or `BIVAR` and passed down through long argument lists rather than local allocation.

## Key conventions

- Keep edits compatible with fixed-form Fortran 77. This codebase uses column-1 `C` comments, `+` continuation lines, and the existing `IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)` convention throughout the numerical routines.
- Keep the hard size limits synchronized in both menu branches. `UNIVAR` and `BIVAR` each declare `PARAMETER(MVAR=4, NDAT=500, IBIN=50)`; if you raise dataset or variable limits, update both places together.
- Input is fixed-width, not free-form. `DATAIN` uses `FORMAT(10(I4,F10.3))` for Kaplan-Meier data and `FORMAT(I4,10(I4,F10.3))` for two-sample data; `DATREG` uses `FORMAT(I4,11F10.3)` for bivariate data.
- If `MVAR` grows past 10, you must widen the `DATAIN` and `DATREG` read formats as described in Appendix A7 of the manual. If data needs more than three decimal places, update the relevant formats and reduce `CONST` in `CENS` so it retains at least two more decimal places than the data.
- `CHARACTER*9` is a real constraint across filenames, variable names, group labels, command files, and output files. Keep those identifiers within 9 characters unless you are intentionally widening the declarations across the input/output code.
- Do not "simplify away" the upper-limit sign normalization in helpers like `AARRAY` and `XVAR`. Those routines flip upper-limit data through `ISIGN` so the downstream statistics can reuse the lower-limit logic.
- Two-sample mode is intentionally bundled: `TWOST` produces the full set of Gehan (permutation and hypergeometric), logrank, Peto-Peto, and Peto-Prentice results rather than asking callers to choose one test at a time.
- When automating runs, prefer command files plus piped top-level menu choices over trying to answer every prompt interactively. The authoritative command-file layouts are in `UNIVAR`, `BIVAR`, and the `R3`-`R6` helper routines, with worked examples mirrored in `manual.txt` and `asurv.etc`.
- `README.md` and `AUTHORS.md` are the active project docs. Preserve the current terminology (`UNIVAR` for Kaplan-Meier / two-sample work, `BIVAR` for correlation / regression work), and note that contributors are expected to discuss substantial changes with the owners before submitting them.
