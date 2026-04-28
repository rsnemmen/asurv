"""Input validation for the ASURV Python API."""

from __future__ import annotations

import pandas as pd

from asurv import AsurvValidationError

_MAX_ROWS = 500
_MAX_VARIABLES = 4
_MAX_BIVAR_INDEPENDENT = 2
_LEGACY_NAME_LIMIT = 9

METHOD_TO_NUMBER: dict[str, int] = {
    "cox": 1,
    "kendall": 2,
    "spearman": 3,
    "em": 4,
    "buckley_james": 5,
    "schmitt": 6,
}

_METHOD_ALIASES: dict[str, str] = {
    "bj": "buckley_james",
    "buckley-james": "buckley_james",
    "buckley_james": "buckley_james",
}


def normalize_method(name: str) -> str:
    normalized = name.strip().lower().replace("-", "_")
    normalized = _METHOD_ALIASES.get(normalized, normalized)
    if normalized not in METHOD_TO_NUMBER:
        raise AsurvValidationError(f"unsupported method: {name!r}")
    return normalized


def legacy_name(value: str, label: str) -> str:
    """Truncate or validate a name to the 9-char Fortran CHARACTER*9 limit."""
    if not value.strip():
        raise AsurvValidationError(f"{label} must not be blank")
    if len(value) > _LEGACY_NAME_LIMIT:
        return value[:_LEGACY_NAME_LIMIT]
    return value


def validate_title(title: str) -> None:
    if not title.strip():
        raise AsurvValidationError("title must not be blank")
    if len(title) > 80:
        raise AsurvValidationError("title must be 80 characters or fewer")


def validate_km_df(
    df: pd.DataFrame, value_col: str, upper_limit_col: str
) -> None:
    if value_col not in df.columns:
        raise AsurvValidationError(f"column {value_col!r} not found in DataFrame")
    if upper_limit_col not in df.columns:
        raise AsurvValidationError(
            f"column {upper_limit_col!r} not found in DataFrame"
        )
    if len(df) == 0:
        raise AsurvValidationError("DataFrame is empty")
    if len(df) > _MAX_ROWS:
        raise AsurvValidationError(
            f"DataFrame has {len(df)} rows; ASURV limit is {_MAX_ROWS}"
        )
    try:
        pd.to_numeric(df[value_col])
    except (ValueError, TypeError) as exc:
        raise AsurvValidationError(
            f"column {value_col!r} contains non-numeric values"
        ) from exc


def validate_twosample_df(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    upper_limit_col: str,
    groups: tuple[object, object],
) -> None:
    for col in (value_col, group_col, upper_limit_col):
        if col not in df.columns:
            raise AsurvValidationError(f"column {col!r} not found in DataFrame")
    if len(df) == 0:
        raise AsurvValidationError("DataFrame is empty")
    if len(df) > _MAX_ROWS:
        raise AsurvValidationError(
            f"DataFrame has {len(df)} rows; ASURV limit is {_MAX_ROWS}"
        )
    if groups[0] == groups[1]:
        raise AsurvValidationError("the two groups must be different")
    data_groups = set(df[group_col].unique())
    for g in groups:
        if g not in data_groups:
            raise AsurvValidationError(
                f"group {g!r} specified in 'groups' does not appear in column {group_col!r}"
            )
    try:
        pd.to_numeric(df[value_col])
    except (ValueError, TypeError) as exc:
        raise AsurvValidationError(
            f"column {value_col!r} contains non-numeric values"
        ) from exc


def validate_bivar_df(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    y_upper_col: str | None,
    x_upper_col: str | None,
) -> None:
    required = [x_col, y_col]
    if y_upper_col is not None:
        required.append(y_upper_col)
    if x_upper_col is not None:
        required.append(x_upper_col)
    for col in required:
        if col not in df.columns:
            raise AsurvValidationError(f"column {col!r} not found in DataFrame")
    if len(df) == 0:
        raise AsurvValidationError("DataFrame is empty")
    if len(df) > _MAX_ROWS:
        raise AsurvValidationError(
            f"DataFrame has {len(df)} rows; ASURV limit is {_MAX_ROWS}"
        )
    for col in (x_col, y_col):
        try:
            pd.to_numeric(df[col])
        except (ValueError, TypeError) as exc:
            raise AsurvValidationError(
                f"column {col!r} contains non-numeric values"
            ) from exc


def validate_bivar_methods(
    methods: list[str], x_upper_col: str | None, n_independent: int
) -> None:
    method_set = set(methods)
    has_x_censoring = x_upper_col is not None
    if n_independent > 1:
        allowed = {"cox", "em", "buckley_james"}
        invalid = method_set - allowed
        if invalid:
            raise AsurvValidationError(
                "multivariate BIVAR only supports cox, em, and buckley-james; "
                f"got: {', '.join(sorted(invalid))}"
            )
        return
    if has_x_censoring:
        allowed = {"kendall", "spearman", "schmitt"}
        invalid = method_set - allowed
        if invalid:
            raise AsurvValidationError(
                "when X is censored, only kendall, spearman, and schmitt are supported; "
                f"got: {', '.join(sorted(invalid))}"
            )
