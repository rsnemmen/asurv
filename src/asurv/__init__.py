"""ASURV: Python interface to the Fortran survival-analysis backend."""

from __future__ import annotations

from asurv.km import KMResult, km
from asurv.twosample import TestStat, TwoSampleResult, twosample
from asurv.bivar import (
    BHKResult,
    BJResult,
    BivariateResult,
    CoxResult,
    EMResult,
    SpearmanResult,
    TwoKMResult,
    bivar,
)

__all__ = [
    # high-level functions
    "km",
    "twosample",
    "bivar",
    # result types
    "KMResult",
    "TwoSampleResult",
    "TestStat",
    "BivariateResult",
    "CoxResult",
    "BHKResult",
    "SpearmanResult",
    "EMResult",
    "BJResult",
    "TwoKMResult",
    # errors
    "AsurvError",
    "AsurvNotFoundError",
    "AsurvValidationError",
    "AsurvExecutionError",
    "AsurvParseError",
]


class AsurvError(Exception):
    """Base exception for all ASURV errors."""


class AsurvNotFoundError(AsurvError):
    """Raised when the asurv executable cannot be located."""


class AsurvValidationError(AsurvError):
    """Raised for invalid inputs before invoking the backend."""


class AsurvExecutionError(AsurvError):
    """Raised when the asurv executable exits with an error."""


class AsurvParseError(AsurvError):
    """Raised when the output file cannot be parsed.

    The raw_output attribute carries the full text for debugging.
    """

    def __init__(self, message: str, raw_output: str = "") -> None:
        super().__init__(message)
        self.raw_output = raw_output
