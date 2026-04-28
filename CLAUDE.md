# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make          # build the asurv executable
make clean    # remove executable
make distclean  # remove executable and compiler scratch files (*.o, *.mod, fort.*)

pip install -e ".[plot,dev]"   # install the Python package (editable, with extras)
```

The Makefile uses `gfortran -std=legacy -fallow-argument-mismatch`. Equivalent manual command:

```bash
gfortran -std=legacy -fallow-argument-mismatch -o asurv asurv.f
```

The Python package depends on the Fortran binary at runtime; `make` must succeed before integration tests will run.

## Tests

Two parallel suites — both must pass:

```bash
# Fortran (authoritative for numerical correctness)
make test              # all smoke + regression tests
make smoke             # output sanity only
make regression        # exact-output diff vs asurv.etc
make {smoke,regression}-gal{1,2,3}  # individual cases

# Python
pytest                            # full suite (53 tests, requires ./asurv built)
pytest -m "not integration"       # offline subset (parsers + validation, no binary needed)
```

The Fortran regression tests extract data/command/expected sections from `asurv.etc` into a temp directory, run the executable via piped stdin, and diff the output. The Python integration tests cover the same three example workflows but assert against parsed result fields (e.g. `result.cox.chi2`) instead of raw text — they catch parser regressions in addition to backend regressions. Numerical correctness is owned by the Fortran suite; treat the Python integration tests as wrapper-correctness checks.

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

**Python package** (`src/asurv/`): the primary user interface. Wraps the Fortran binary in three layers:

1. **Input translation** (`_ind.py`, `_validate.py`, `_commands.py`) — pandas DataFrame + boolean censoring columns → ASURV IND integer codes + fixed-width data text + command-file lines
2. **Backend invocation** (`_backend.py`) — temp workspace, executable discovery (`ASURV_BIN` env var → `./asurv` in cwd → `PATH`), `subprocess.run` in batch mode
3. **Output parsing** (`_parsers.py`) — section-based regex on `.out` text → frozen dataclasses (`KMResult`, `TwoSampleResult`, `BivariateResult`, etc.)

Public API in `src/asurv/__init__.py`: `asurv.km()`, `asurv.twosample()`, `asurv.bivar()`. Result dataclasses carry `.raw_output` (full text), `.summary()` (human digest), and `.plot()` (lazy matplotlib).

CLI entry point: `asurv-py {km,twost,bivar} ...` (defined in `src/asurv/cli.py`). The legacy `asurv_cli.py` at the repo root is a 3-line backward-compatibility shim — do not add logic to it.

**Wrap, don't rewrite.** Do not extend `asurv.f` to add Python-facing features (structured output, JSON, etc.). The Fortran is the source of numerical truth and must remain bit-for-bit reproducible against `asurv.etc`. New features belong in the Python layer.

**When the Fortran output format changes:** update parsers in `src/asurv/_parsers.py` only. Anchor strings to know about: `# OF DATA POINTS`, `VARIABLE RANGE      KM ESTIMATOR   ERROR`, `PERCENTILES`, `MEAN=`, `DIFFERENTIAL KM ESTIMATOR`, `TEST STATISTIC`, `PROBABILITY`, `GLOBAL CHI SQUARE`, `INTERCEPT COEFF`, `SLOPE COEFF`. Test fixtures are extracted from `asurv.etc` by `tests/conftest.py`.

**Automating runs without the Python wrapper:** prefer command files plus piped top-level menu choices over scripting every interactive prompt. The authoritative command-file formats are defined in `UNIVAR`, `BIVAR`, and `R3`–`R6` in `asurv.f`; worked examples are in `manual.txt` and `asurv.etc`.

## Key Fortran Conventions

- Fixed-form Fortran 77: column-1 `C` comments, `+` continuation lines, `IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)` throughout numerical routines.
- Hard size limits declared as `PARAMETER(MVAR=4, NDAT=500, IBIN=50)` in **both** `UNIVAR` and `BIVAR`. If raising limits, update both.
- Input is fixed-width: `DATAIN` uses `FORMAT(10(I4,F10.3))` for KM data and `FORMAT(I4,10(I4,F10.3))` for two-sample data; `DATREG` uses `FORMAT(I4,11F10.3)` for bivariate data. If `MVAR` grows past 10, widen these formats per Appendix A7 of the manual. If data needs more than three decimal places, update the relevant formats and reduce `CONST` in `CENS` so it retains at least two more decimal places than the data.
- `CHARACTER*9` is a hard limit for filenames, variable names, group labels, command files, and output files.
- Do not remove the upper-limit sign normalization in `AARRAY` and `XVAR` — these routines flip upper-limit data via `ISIGN` so downstream statistics can reuse lower-limit logic.
- Work arrays are allocated in `UNIVAR`/`BIVAR` and passed through long argument lists rather than local allocation.
- `TWOST` always computes all tests (Gehan permutation, Gehan hypergeometric, logrank, Peto-Peto, Peto-Prentice) — this is intentional, not redundant.

### Quirks the Python wrapper hides (re-introduce only if bypassing the wrapper)

- **Output files use `STATUS='NEW'`** — the named `.out` file must not pre-exist, or the open fails. The Python `Workspace` context manager creates a fresh temp dir per call.
- **Schmitt bootstrap requires a *negative* random seed.** `bivar_command_lines` in `_commands.py` silently negates positive values.
- **`OUTPUT='         '` (9 spaces) routes results to stdout** instead of a file — mostly used by the interactive path.
- **Filenames truncated to 9 chars.** The Python wrapper always copies user data to `data.dat` in the workspace and uses short names (`km.out`, `twost.out`, `bivar.out`) so user paths never hit the limit.
- **Command-file lines parsed by `DATA2`** are 4-wide fields (`4A1` / `12A1` / `20A1`); integers must fit in 4 columns.

## Sample Data

`asurv.etc` contains embedded sections delimited by `****` headers:
- `gal1.dat` / `gal1.com` / `gal1.out` — Kaplan-Meier example (IR Luminosities)
- `gal2.dat` / `gal2.com` / `gal2.out` — two-sample test example
- `gal3.dat` / `gal3.com` / `gal3.out` — bivariate correlation/regression example

Both test suites consume these as fixtures: the Makefile uses awk pipelines (`Makefile:51-118`); the Python suite uses a regex split in `tests/conftest.py`.

## Repo layout

```
asurv.f, Makefile, asurv.etc, manual.txt   Fortran sources + sample data
pyproject.toml                              Python package config
src/asurv/                                  Python package (wrappers, parsers, results)
tests/                                      pytest suite
asurv_cli.py                                3-line backward-compat shim → src/asurv/cli.py
```
