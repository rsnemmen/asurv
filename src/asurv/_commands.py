"""Generate ASURV command-file lines and fixed-width data file text.

Command-file formats are derived from UVCMD (asurv.f:9528) and BVCMD
(asurv.f:2304).  The fixed-width data formats come from DATAIN (asurv.f:3439)
and DATREG (asurv.f:2991).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from asurv._validate import METHOD_TO_NUMBER, legacy_name


# ---------------------------------------------------------------------------
# Data file serialisers
# ---------------------------------------------------------------------------

def km_data_text(df: pd.DataFrame, value_col: str, ind: np.ndarray) -> str:
    """Serialise KM data to FORMAT(10(I4,F10.3)) fixed-width text."""
    lines = []
    for ind_val, x in zip(ind, df[value_col]):
        lines.append(f"{int(ind_val):4d}{float(x):10.3f}")
    return "\n".join(lines) + "\n"


def twost_data_text(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    ind: np.ndarray,
    group_map: dict[object, int],
) -> str:
    """Serialise two-sample data to FORMAT(I4,10(I4,F10.3)) fixed-width text."""
    lines = []
    for group_val, ind_val, x in zip(df[group_col], ind, df[value_col]):
        g = group_map[group_val]
        lines.append(f"{g:4d}{int(ind_val):4d}{float(x):10.3f}")
    return "\n".join(lines) + "\n"


def bivar_data_text(
    df: pd.DataFrame,
    x_cols: list[str],
    y_col: str,
    ind: np.ndarray,
) -> str:
    """Serialise bivariate data to FORMAT(I4,11F10.3) fixed-width text."""
    lines = []
    for i, row in enumerate(df.itertuples(index=False)):
        x_vals = "".join(f"{float(getattr(row, c)):10.3f}" for c in x_cols)
        y_val = float(getattr(row, y_col))
        lines.append(f"{int(ind[i]):4d}{x_vals}{y_val:10.3f}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Command-file line generators
# ---------------------------------------------------------------------------

def km_command_lines(
    data_filename: str,
    title: str,
    nvar: int,
    column: int,
    variable_name: str,
    output_filename: str,
    full_km: bool = True,
    differential: tuple[float, int, float] | None = None,
    print_data: bool = False,
) -> list[str]:
    """Return the ordered list of lines for a KM command file."""
    vname = legacy_name(variable_name, "variable name")
    kdiff = 1 if differential else 0
    diff_start = differential[0] if differential else 0.0
    diff_bins = differential[1] if differential else 1
    diff_size = differential[2] if differential else 1.0
    return [
        data_filename,
        title,
        str(nvar),
        str(column),
        vname,
        "1" if full_km else "0",
        str(kdiff),
        f"{diff_start:.3f}",
        str(diff_bins),
        f"{diff_size:.3f}",
        "1" if print_data else "0",
        output_filename,
    ]


def twost_command_lines(
    data_filename: str,
    title: str,
    nvar: int,
    column: int,
    variable_name: str,
    group_ids: list[int],
    first_group_id: int,
    second_group_id: int,
    first_group_name: str,
    second_group_name: str,
    output_filename: str,
    details: bool = False,
    full_km: bool = False,
    differential: tuple[float, int, float] | None = None,
    print_data: bool = False,
) -> list[str]:
    """Return the ordered list of lines for a TWOST command file."""
    vname = legacy_name(variable_name, "variable name")
    g1name = legacy_name(first_group_name, "first group name")
    g2name = legacy_name(second_group_name, "second group name")
    kdiff = 1 if differential else 0
    diff_start = differential[0] if differential else 0.0
    diff_bins = differential[1] if differential else 1
    diff_size = differential[2] if differential else 1.0

    lines = [
        data_filename,
        title,
        str(nvar),
        str(column),
        vname,
        str(len(group_ids)),
        "   ".join(str(g) for g in group_ids),
        "   ".join([
            str(first_group_id),
            str(second_group_id),
            "1" if details else "0",
            "1" if full_km else "0",
        ]),
        str(kdiff),
    ]
    if kdiff:
        lines += [f"{diff_start:.3f}", str(diff_bins), f"{diff_size:.3f}"]
    lines += [
        "1" if print_data else "0",
        g1name,
        g2name,
        output_filename,
    ]
    return lines


def _render_double_a9(first: str, second: str) -> str:
    return f"{first:<9}{second if len(second) == 9 else second.rjust(9)}"


def bivar_command_lines(
    data_filename: str,
    title: str,
    n_independent: int,
    column: int,
    x_name: str,
    y_name: str,
    methods: list[str],
    output_filename: str,
    print_data: bool = False,
    spearman_print_ranks: bool = False,
    em_tolerance: float = 1e-5,
    em_max_iterations: int = 50,
    em_initial_coefficients: list[float] | None = None,
    bj_tolerance: float = 1e-5,
    bj_max_iterations: int = 50,
    schmitt_x_bins: int | None = None,
    schmitt_y_bins: int | None = None,
    schmitt_use_bin_geometry: bool = False,
    schmitt_tolerance: float = 1e-5,
    schmitt_max_iterations: int = 50,
    schmitt_x_bin_size: float | None = None,
    schmitt_y_bin_size: float | None = None,
    schmitt_x_origin: float | None = None,
    schmitt_y_origin: float | None = None,
    schmitt_print_2d_km: bool = False,
    schmitt_bootstrap_iterations: int = 0,
    schmitt_random_seed: int | None = None,
) -> list[str]:
    """Return the ordered list of lines for a BIVAR command file."""
    from asurv import AsurvValidationError

    xname = legacy_name(x_name, "x name")
    yname = legacy_name(y_name, "y name")

    lines = [
        title,
        data_filename,
        "   ".join([str(n_independent), str(column), str(len(methods))]),
        "   ".join(str(METHOD_TO_NUMBER[m]) for m in methods),
        _render_double_a9(xname, yname),
        "1" if print_data else "0",
        output_filename,
    ]

    for method in methods:
        if method == "spearman":
            lines.append("1" if spearman_print_ranks else "0")
        elif method == "em":
            expected = n_independent + 2
            coeffs = em_initial_coefficients or [0.0] * expected
            if len(coeffs) != expected:
                raise AsurvValidationError(
                    f"EM needs {expected} initial coefficients, got {len(coeffs)}"
                )
            lines.append(f"{em_tolerance:.3E}")
            lines.append("".join(f"{v:10.3f}" for v in coeffs))
            lines.append(str(em_max_iterations))
        elif method == "buckley_james":
            lines.append(f"{bj_tolerance:.3E}")
            lines.append(str(bj_max_iterations))
        elif method == "schmitt":
            if schmitt_x_bins is None or schmitt_y_bins is None:
                raise AsurvValidationError(
                    "schmitt method requires schmitt_x_bins and schmitt_y_bins"
                )
            lines.append(f"{schmitt_x_bins:4d}{schmitt_y_bins:4d}")
            lines.append("1" if schmitt_use_bin_geometry else "0")
            lines.append(f"{schmitt_tolerance:.3E}")
            lines.append(str(schmitt_max_iterations))
            if schmitt_use_bin_geometry:
                if None in (schmitt_x_bin_size, schmitt_y_bin_size,
                            schmitt_x_origin, schmitt_y_origin):
                    raise AsurvValidationError(
                        "schmitt bin geometry requires x/y bin sizes and origins"
                    )
                lines.append(
                    f"{schmitt_x_bin_size:10.3f}{schmitt_y_bin_size:10.3f}"
                )
                lines.append(
                    f"{schmitt_x_origin:10.3f}{schmitt_y_origin:10.3f}"
                )
            else:
                lines.append(f"{0.0:10.3f}{0.0:10.3f}")
                lines.append(f"{0.0:10.3f}{0.0:10.3f}")
            lines.append("1" if schmitt_print_2d_km else "0")
            if schmitt_bootstrap_iterations < 0:
                raise AsurvValidationError(
                    "schmitt_bootstrap_iterations must be non-negative"
                )
            lines.append(str(schmitt_bootstrap_iterations))
            if schmitt_bootstrap_iterations > 0:
                if schmitt_random_seed is None:
                    raise AsurvValidationError(
                        "schmitt bootstrap requires schmitt_random_seed"
                    )
                # Fortran requires a negative seed; negate silently if positive
                seed = schmitt_random_seed
                if seed > 0:
                    seed = -seed
                lines.append(str(seed))

    return lines
