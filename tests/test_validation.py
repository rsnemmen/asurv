"""Unit tests for _validate.py — all tests are offline (no Fortran binary)."""

from __future__ import annotations

import pandas as pd
import pytest

from asurv import AsurvValidationError
from asurv._validate import (
    normalize_method,
    validate_bivar_df,
    validate_km_df,
    validate_title,
    validate_twosample_df,
)


class TestTitle:
    def test_blank(self) -> None:
        with pytest.raises(AsurvValidationError, match="blank"):
            validate_title("")

    def test_too_long(self) -> None:
        with pytest.raises(AsurvValidationError, match="80"):
            validate_title("x" * 81)

    def test_valid(self) -> None:
        validate_title("A valid title")


class TestNormalizeMethod:
    def test_aliases(self) -> None:
        assert normalize_method("bj") == "buckley_james"
        assert normalize_method("buckley-james") == "buckley_james"

    def test_invalid(self) -> None:
        with pytest.raises(AsurvValidationError, match="unsupported"):
            normalize_method("bogus")

    def test_cox(self) -> None:
        assert normalize_method("cox") == "cox"


class TestValidateKmDf:
    def _good(self) -> pd.DataFrame:
        return pd.DataFrame({"v": [1.0, 2.0], "ul": [False, True]})

    def test_missing_value_col(self) -> None:
        df = self._good().rename(columns={"v": "x"})
        with pytest.raises(AsurvValidationError, match="'v'"):
            validate_km_df(df, "v", "ul")

    def test_missing_upper_limit_col(self) -> None:
        df = self._good()
        with pytest.raises(AsurvValidationError, match="'flag'"):
            validate_km_df(df, "v", "flag")

    def test_empty(self) -> None:
        with pytest.raises(AsurvValidationError, match="empty"):
            validate_km_df(pd.DataFrame({"v": [], "ul": []}), "v", "ul")

    def test_too_many_rows(self) -> None:
        df = pd.DataFrame({"v": [1.0] * 501, "ul": [False] * 501})
        with pytest.raises(AsurvValidationError, match="500"):
            validate_km_df(df, "v", "ul")

    def test_valid(self) -> None:
        validate_km_df(self._good(), "v", "ul")


class TestValidateTwoSampleDf:
    def _good(self) -> pd.DataFrame:
        return pd.DataFrame({
            "value": [1.0, 2.0, 3.0, 4.0],
            "group": ["A", "A", "B", "B"],
            "ul": [False, True, False, False],
        })

    def test_same_groups(self) -> None:
        df = self._good()
        with pytest.raises(AsurvValidationError, match="different"):
            validate_twosample_df(df, "value", "group", "ul", ("A", "A"))

    def test_group_not_in_data(self) -> None:
        df = self._good()
        with pytest.raises(AsurvValidationError, match="'C'"):
            validate_twosample_df(df, "value", "group", "ul", ("A", "C"))

    def test_valid(self) -> None:
        df = self._good()
        validate_twosample_df(df, "value", "group", "ul", ("A", "B"))


class TestValidateBivarDf:
    def _good(self) -> pd.DataFrame:
        return pd.DataFrame({
            "x": [1.0, 2.0],
            "y": [2.0, 3.0],
            "yu": [False, True],
        })

    def test_missing_x(self) -> None:
        df = self._good()
        with pytest.raises(AsurvValidationError, match="'z'"):
            validate_bivar_df(df, "z", "y", "yu", None)

    def test_missing_y_upper(self) -> None:
        df = self._good()
        with pytest.raises(AsurvValidationError, match="'noexist'"):
            validate_bivar_df(df, "x", "y", "noexist", None)

    def test_valid(self) -> None:
        df = self._good()
        validate_bivar_df(df, "x", "y", "yu", None)
