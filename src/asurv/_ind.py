"""Translate user-facing boolean censoring flags into ASURV IND integer codes.

ASURV censoring conventions
----------------------------
KM / TWOST (DATAIN):
  IND > 0  detected (lower limit / normal measurement)
  IND < 0  upper limit (censored from above)
  The Fortran uses the sign: detected = positive, upper limit = negative.
  In practice the code checks ``ISIGN`` and uses magnitude, so any positive
  integer means "detected" and any negative integer means "upper limit".
  We use 0 = detected, -1 = upper limit to match the sample files.

BIVAR (DATREG):
  0   both detected
  1   Y lower limit  (not applicable here — we only expose upper limits)
 -1   Y upper limit
  2   X lower limit
 -2   X upper limit
  3   both lower limits
 -3   Y-upper / X-lower  (treated as both upper here — simplification)
  4   Y lower / X upper  (not exposed)
 -4   both upper limits (our convention when both are upper limits)

For the common case (upper limits only) astronomers will pass a single boolean
column.  For bivariate data with censoring in X and/or Y the caller passes
separate boolean columns.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def km_ind(df: pd.DataFrame, upper_limit_col: str) -> np.ndarray:
    """Return a 1-D integer array of KM IND codes (0 = detected, -1 = upper limit)."""
    flags = df[upper_limit_col].astype(bool).to_numpy()
    return np.where(flags, -1, 0).astype(int)


def twost_ind(df: pd.DataFrame, upper_limit_col: str) -> np.ndarray:
    """Same as km_ind but for two-sample data."""
    return km_ind(df, upper_limit_col)


def bivar_ind(
    df: pd.DataFrame,
    y_upper_col: str | None = None,
    x_upper_col: str | None = None,
) -> np.ndarray:
    """Return a 1-D integer array of BIVAR IND codes.

    IND = 0   both detected
    IND = -1  Y upper limit only
    IND = -2  X upper limit only
    IND = -4  both upper limits
    """
    n = len(df)
    y_upper = (
        df[y_upper_col].astype(bool).to_numpy()
        if y_upper_col is not None
        else np.zeros(n, dtype=bool)
    )
    x_upper = (
        df[x_upper_col].astype(bool).to_numpy()
        if x_upper_col is not None
        else np.zeros(n, dtype=bool)
    )
    codes = np.zeros(n, dtype=int)
    codes[y_upper & ~x_upper] = -1
    codes[~y_upper & x_upper] = -2
    codes[y_upper & x_upper] = -4
    return codes
