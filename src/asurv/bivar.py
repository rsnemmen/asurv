"""Bivariate correlation and regression for censored data."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class CoxResult:
    chi2: float
    dof: int
    probability: float

    def __repr__(self) -> str:
        return f"CoxResult(chi2={self.chi2:.4f}, dof={self.dof}, p={self.probability:.4f})"


@dataclass(frozen=True)
class BHKResult:
    z: float
    probability: float

    def __repr__(self) -> str:
        return f"BHKResult(z={self.z:.4f}, p={self.probability:.4f})"


@dataclass(frozen=True)
class SpearmanResult:
    rho: float
    z: float
    probability: float

    def __repr__(self) -> str:
        return f"SpearmanResult(rho={self.rho:.4f}, z={self.z:.4f}, p={self.probability:.4f})"


@dataclass(frozen=True)
class EMResult:
    intercept: float
    intercept_err: float
    slopes: list[tuple[float, float]]
    sigma: float
    iterations: int

    def __repr__(self) -> str:
        s = self.slopes[0] if self.slopes else ("?", "?")
        return (
            f"EMResult(intercept={self.intercept:.4f}±{self.intercept_err:.4f}, "
            f"slope={s[0]:.4f}±{s[1]:.4f}, sigma={self.sigma:.4f})"
        )


@dataclass(frozen=True)
class BJResult:
    intercept: float
    intercept_err: float
    slope: float
    slope_err: float

    def __repr__(self) -> str:
        return (
            f"BJResult(intercept={self.intercept:.4f}±{self.intercept_err:.4f}, "
            f"slope={self.slope:.4f}±{self.slope_err:.4f})"
        )


@dataclass(frozen=True)
class TwoKMResult:
    slope: float
    intercept: float

    def __repr__(self) -> str:
        return f"TwoKMResult(slope={self.slope:.4f}, intercept={self.intercept:.4f})"


@dataclass(frozen=True)
class BivariateResult:
    """Results of ASURV bivariate correlation and regression methods."""

    n: int
    upper_limits: dict[str, int]
    cox: CoxResult | None
    bhk: BHKResult | None
    spearman: SpearmanResult | None
    em: EMResult | None
    bj: BJResult | None
    twokm: TwoKMResult | None
    raw_output: str

    def summary(self) -> str:
        lines = [
            "Bivariate Correlation & Regression",
            f"  N = {self.n}",
            f"  Upper limits — Y: {self.upper_limits.get('y', 0)}  "
            f"X: {self.upper_limits.get('x', 0)}  "
            f"Both: {self.upper_limits.get('both', 0)}",
        ]
        if self.cox:
            lines.append(
                f"  Cox:      χ²={self.cox.chi2:.4f}  dof={self.cox.dof}  p={self.cox.probability:.4f}"
            )
        if self.bhk:
            lines.append(
                f"  BHK:      Z={self.bhk.z:.4f}  p={self.bhk.probability:.4f}"
            )
        if self.spearman:
            lines.append(
                f"  Spearman: ρ={self.spearman.rho:.4f}  Z={self.spearman.z:.4f}  p={self.spearman.probability:.4f}"
            )
        if self.em:
            slope_str = "  ".join(f"{s:.4f}±{e:.4f}" for s, e in self.em.slopes)
            lines.append(
                f"  EM:       intercept={self.em.intercept:.4f}±{self.em.intercept_err:.4f}  "
                f"slope(s): {slope_str}  σ={self.em.sigma:.4f}  iter={self.em.iterations}"
            )
        if self.bj:
            lines.append(
                f"  B-J:      intercept={self.bj.intercept:.4f}±{self.bj.intercept_err:.4f}  "
                f"slope={self.bj.slope:.4f}±{self.bj.slope_err:.4f}"
            )
        if self.twokm:
            lines.append(
                f"  TWOKM:    slope={self.twokm.slope:.4f}  intercept={self.twokm.intercept:.4f}"
            )
        return "\n".join(lines)

    def plot(
        self,
        df: pd.DataFrame | None = None,
        x: str | None = None,
        y: str | None = None,
        x_upper: str | None = None,
        y_upper: str | None = None,
        ax: Any = None,
    ) -> Any:
        """Scatter plot of data with regression line overlays.

        Requires matplotlib (``pip install asurv[plot]``).

        Parameters
        ----------
        df:
            The same DataFrame passed to :func:`bivar`.  If omitted, only the
            regression lines are drawn (requires ``x`` limits from the caller).
        x, y:
            Column names for the independent and dependent variables.
        x_upper, y_upper:
            Boolean column names for X and Y censoring (``True`` = upper limit).
        ax:
            Existing matplotlib Axes.
        """
        from asurv._plots import plot_scatter_regression

        x_vals = df[x].tolist() if df is not None and x else None
        y_vals = df[y].tolist() if df is not None and y else None
        xu = df[x_upper].tolist() if df is not None and x_upper else None
        yu = df[y_upper].tolist() if df is not None and y_upper else None

        em_intercept = self.em.intercept if self.em else None
        em_slopes = [s for s, _ in self.em.slopes] if self.em else None
        bj_int = self.bj.intercept if self.bj else None
        bj_sl = self.bj.slope if self.bj else None

        return plot_scatter_regression(
            x=x_vals or [],
            y=y_vals or [],
            x_upper=xu,
            y_upper=yu,
            intercept=em_intercept,
            slopes=em_slopes,
            bj_intercept=bj_int,
            bj_slope=bj_sl,
            ax=ax,
            x_label=x or "X",
            y_label=y or "Y",
        )

    def __repr__(self) -> str:
        present = [
            name
            for name in ("cox", "bhk", "spearman", "em", "bj", "twokm")
            if getattr(self, name) is not None
        ]
        return f"BivariateResult(n={self.n}, methods={present})"


def bivar(
    df: pd.DataFrame,
    x: str,
    y: str,
    methods: list[str] | str,
    y_upper: str | None = None,
    x_upper: str | None = None,
    title: str = "Bivariate Analysis",
    x_name: str | None = None,
    y_name: str | None = None,
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
) -> BivariateResult:
    """Run ASURV bivariate correlation and/or regression on a censored dataset.

    Parameters
    ----------
    df:
        DataFrame containing the data.
    x:
        Column name of the independent variable.
    y:
        Column name of the dependent variable.
    methods:
        One or more of ``"cox"``, ``"kendall"``, ``"spearman"``, ``"em"``,
        ``"buckley-james"`` (alias ``"bj"``), ``"schmitt"``.
    y_upper:
        Boolean column name: ``True`` = Y is an upper limit.
    x_upper:
        Boolean column name: ``True`` = X is an upper limit.
    title:
        Analysis title (≤80 characters).
    x_name:
        Label for X (≤9 characters; defaults to column name).
    y_name:
        Label for Y (≤9 characters; defaults to column name).

    Returns
    -------
    BivariateResult
    """
    from asurv._backend import run
    from asurv._commands import bivar_command_lines, bivar_data_text
    from asurv._ind import bivar_ind
    from asurv._parsers import parse_bivar
    from asurv._validate import (
        normalize_method,
        validate_bivar_df,
        validate_bivar_methods,
        validate_title,
    )

    if isinstance(methods, str):
        methods = [methods]

    validate_title(title)
    validate_bivar_df(df, x, y, y_upper, x_upper)

    norm_methods = [normalize_method(m) for m in methods]
    if len(set(norm_methods)) != len(norm_methods):
        from asurv import AsurvValidationError
        raise AsurvValidationError("duplicate methods are not allowed")
    validate_bivar_methods(norm_methods, x_upper, n_independent=1)

    xname = (x_name or x)[:9]
    yname = (y_name or y)[:9]

    ind = bivar_ind(df, y_upper_col=y_upper, x_upper_col=x_upper)
    data_text = bivar_data_text(df, [x], y, ind)

    cmd = bivar_command_lines(
        data_filename="data.dat",
        title=title,
        n_independent=1,
        column=1,
        x_name=xname,
        y_name=yname,
        methods=norm_methods,
        output_filename="bivar.out",
        print_data=print_data,
        spearman_print_ranks=spearman_print_ranks,
        em_tolerance=em_tolerance,
        em_max_iterations=em_max_iterations,
        em_initial_coefficients=em_initial_coefficients,
        bj_tolerance=bj_tolerance,
        bj_max_iterations=bj_max_iterations,
        schmitt_x_bins=schmitt_x_bins,
        schmitt_y_bins=schmitt_y_bins,
        schmitt_use_bin_geometry=schmitt_use_bin_geometry,
        schmitt_tolerance=schmitt_tolerance,
        schmitt_max_iterations=schmitt_max_iterations,
        schmitt_x_bin_size=schmitt_x_bin_size,
        schmitt_y_bin_size=schmitt_y_bin_size,
        schmitt_x_origin=schmitt_x_origin,
        schmitt_y_origin=schmitt_y_origin,
        schmitt_print_2d_km=schmitt_print_2d_km,
        schmitt_bootstrap_iterations=schmitt_bootstrap_iterations,
        schmitt_random_seed=schmitt_random_seed,
    )
    raw = run("BIVAR", data_text, cmd,
              command_filename="bivar.com", output_filename="bivar.out")
    return parse_bivar(raw)
