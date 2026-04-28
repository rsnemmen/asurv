# ASURV Copilot Instructions

## Build and verification

Build the Fortran backend, then install the Python package:

```bash
make                            # build ./asurv (gfortran -std=legacy -fallow-argument-mismatch)
pip install -e ".[plot,dev]"    # install the Python wrapper (editable, with extras)
```

Two parallel test suites — both must pass:

```bash
# Fortran (authoritative for numerical correctness)
make test                        # all smoke + regression tests
make {smoke,regression}-gal{1,2,3}  # individual cases

# Python
pytest                           # 53 tests (requires ./asurv built)
pytest -m "not integration"      # offline subset; no binary needed
```

The Fortran regression tests extract `gal1`/`gal2`/`gal3` data from `asurv.etc`, run the binary via piped stdin, and diff the output. The Python integration tests assert against parsed result fields (e.g. `result.cox.chi2`) and catch parser regressions in addition to backend regressions. Numerical correctness is owned by the Fortran suite.

## High-level architecture

`asurv.f` is the whole Fortran application: one fixed-form Fortran 77 source file containing the main menu plus all statistical subroutines. The old multi-file build list is preserved only as historical reference inside `asurv.etc`.

**Fortran dispatch:**
- `MAIN` → `UNIVAR` (menu choice 1) or `BIVAR` (menu choice 2).
- `UNIVAR` handles Kaplan-Meier and two-sample tests. Reads data via `DATAIN`, calls `KMESTM` and `TWOST`.
- `BIVAR` handles correlation/regression. Reads data via `DATREG`, dispatches to `COXREG`, `BHK`, `SPRMAN`, `EM`, `BJ`, or `TWOKM`. Method-specific prompts are in `R3`–`R6`.
- Most algorithm routines work on fixed-size work arrays allocated in `UNIVAR`/`BIVAR` and passed through long argument lists.

**Batch modes** (bypass the interactive menu):

```bash
./asurv KM km.com
./asurv TWOST ts.com
./asurv BIVAR bv.com
```

**Python package** (`src/asurv/`): the primary user interface. Wraps the Fortran binary in three layers:

1. **Input translation** (`_ind.py`, `_validate.py`, `_commands.py`) — pandas DataFrame + boolean censoring columns → ASURV IND codes + fixed-width data text + command-file lines
2. **Backend invocation** (`_backend.py`) — temp workspace, executable discovery (`ASURV_BIN` env var → `./asurv` in cwd → PATH), `subprocess.run` in batch mode
3. **Output parsing** (`_parsers.py`) — section-based regex on `.out` text → frozen dataclasses (`KMResult`, `TwoSampleResult`, `BivariateResult`, etc.)

Public API: `asurv.km()`, `asurv.twosample()`, `asurv.bivar()`. CLI entry point: `asurv-py`. The legacy `asurv_cli.py` at the repo root is a 3-line backward-compat shim — do not add logic to it.

**Wrap, don't rewrite.** Do not extend `asurv.f` to add Python-facing features. The Fortran is the source of numerical truth and must remain bit-for-bit reproducible against `asurv.etc`. New features belong in the Python layer.

**When Fortran output format changes:** update `src/asurv/_parsers.py` only. Key anchor strings: `# OF DATA POINTS`, `VARIABLE RANGE      KM ESTIMATOR   ERROR`, `MEAN=`, `TEST STATISTIC`, `PROBABILITY`, `GLOBAL CHI SQUARE`, `INTERCEPT COEFF`, `SLOPE COEFF`.

## Key conventions

- Keep edits compatible with fixed-form Fortran 77: column-1 `C` comments, `+` continuation lines, `IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)` throughout numerical routines.
- Keep the hard size limits synchronized in both menu branches. `UNIVAR` and `BIVAR` each declare `PARAMETER(MVAR=4, NDAT=500, IBIN=50)`; if you raise limits, update both.
- Input is fixed-width: `DATAIN` uses `FORMAT(10(I4,F10.3))` for KM data and `FORMAT(I4,10(I4,F10.3))` for two-sample data; `DATREG` uses `FORMAT(I4,11F10.3)` for bivariate data.
- If `MVAR` grows past 10, widen the `DATAIN`/`DATREG` formats per Appendix A7 of the manual. If data needs more than three decimal places, update formats and reduce `CONST` in `CENS`.
- `CHARACTER*9` is a hard limit for filenames, variable names, group labels, command files, and output files.
- Do not remove the upper-limit sign normalization in `AARRAY` and `XVAR` — those routines flip upper-limit data via `ISIGN` so downstream statistics can reuse lower-limit logic.
- `TWOST` always computes all five tests (Gehan permutation, Gehan hypergeometric, logrank, Peto-Peto, Peto-Prentice) — this is intentional, not redundant.

### Quirks the Python wrapper hides (relevant if bypassing the wrapper)

- **Output files use `STATUS='NEW'`** — the named `.out` file must not pre-exist or the open fails.
- **Schmitt bootstrap requires a negative random seed** — `bivar_command_lines` in `_commands.py` silently negates positive values.
- **`OUTPUT='         '` (9 spaces) routes results to stdout** instead of a file.
- **Filenames are truncated to 9 chars** — the Python wrapper always copies user data to `data.dat` in a fresh temp workspace.
- **Command-file integers are 4-wide** — fields read by `DATA2` use `4A1`/`12A1`/`20A1` formats.

## Repo layout

```
asurv.f, Makefile, asurv.etc, manual.txt   Fortran sources + sample data
pyproject.toml                              Python package config
src/asurv/                                  Python package (wrappers, parsers, results)
tests/                                      pytest suite (conftest.py extracts gal1/2/3 from asurv.etc)
asurv_cli.py                                3-line backward-compat shim → src/asurv/cli.py
```

`README.md` is the active user-facing doc. Preserve existing terminology (`UNIVAR` for KM/two-sample, `BIVAR` for correlation/regression). Discuss substantial changes with the owners before submitting.
