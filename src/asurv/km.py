"""Kaplan-Meier estimator API."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class KMResult:
    """Results of the Kaplan-Meier estimator."""

    n: int
    n_upper_limits: int
    table: pd.DataFrame
    censored_points: list[float]
    percentiles: dict[int, float] | None
    mean: float
    mean_err: float
    differential: pd.DataFrame | None
    raw_output: str

    def summary(self) -> str:
        lines = [
            f"Kaplan-Meier Estimator",
            f"  N = {self.n}  (upper limits: {self.n_upper_limits})",
            f"  Mean = {self.mean:.4f} ± {self.mean_err:.4f}",
        ]
        if self.percentiles:
            pct = "  Percentiles: " + "  ".join(
                f"{k}th={v:.4f}" for k, v in sorted(self.percentiles.items())
            )
            lines.append(pct)
        return "\n".join(lines)

    def plot(self, ax: Any = None, label: str | None = None, **kwargs: Any) -> Any:
        """Plot the survival step function with censoring tick marks.

        Requires matplotlib (``pip install asurv[plot]``).
        """
        from asurv._plots import plot_km

        return plot_km(self.table, self.censored_points, label=label, ax=ax, **kwargs)

    def __repr__(self) -> str:
        return (
            f"KMResult(n={self.n}, n_upper_limits={self.n_upper_limits}, "
            f"mean={self.mean:.4f}±{self.mean_err:.4f})"
        )


def km(
    df: pd.DataFrame,
    value: str,
    upper_limit: str,
    title: str = "KM Estimator",
    variable_name: str | None = None,
    full_km: bool = True,
    differential: tuple[float, int, float] | None = None,
    print_data: bool = False,
) -> KMResult:
    """Run the Kaplan-Meier estimator on a censored dataset.

    Parameters
    ----------
    df:
        DataFrame with at least two columns: the measured values and a boolean
        censoring flag.
    value:
        Name of the column containing measured values.
    upper_limit:
        Name of the boolean column where ``True`` means "upper limit"
        (censored from above).
    title:
        Title string for the analysis (≤80 characters).
    variable_name:
        Label for the variable (≤9 characters; defaults to the column name).
    full_km:
        Whether to print the full KM table in the output (default ``True``).
    differential:
        If not ``None``, compute the differential KM histogram.  Pass a tuple
        ``(start, n_bins, bin_size)``.
    print_data:
        Whether to include the input data in the output.

    Returns
    -------
    KMResult
    """
    from asurv._backend import run
    from asurv._commands import km_command_lines, km_data_text
    from asurv._ind import km_ind
    from asurv._parsers import parse_km
    from asurv._validate import validate_km_df, validate_title

    validate_title(title)
    validate_km_df(df, value, upper_limit)

    vname = (variable_name or value)[:9]
    ind = km_ind(df, upper_limit)
    data_text = km_data_text(df, value, ind)
    cmd = km_command_lines(
        data_filename="data.dat",
        title=title,
        nvar=1,
        column=1,
        variable_name=vname,
        output_filename="km.out",
        full_km=full_km,
        differential=differential,
        print_data=print_data,
    )
    raw = run("KM", data_text, cmd, output_filename="km.out")
    return parse_km(raw)
