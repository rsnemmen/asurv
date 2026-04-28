"""Two-sample tests for censored data."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class TestStat:
    """A single test statistic and its p-value."""

    statistic: float
    probability: float

    def __repr__(self) -> str:
        return f"TestStat(statistic={self.statistic:.4f}, probability={self.probability:.4f})"


@dataclass(frozen=True)
class TwoSampleResult:
    """Results of the ASURV two-sample tests."""

    n: int
    n_upper_limits: int
    n_group1: int
    n_group2: int
    tests: dict[str, TestStat]
    km_group1: "KMResult | None"
    km_group2: "KMResult | None"
    raw_output: str

    def summary(self) -> str:
        lines = [
            "Two-Sample Tests",
            f"  N = {self.n}  (upper limits: {self.n_upper_limits})",
            f"  Group I: {self.n_group1}   Group II: {self.n_group2}",
            "",
            f"  {'Test':<40}  {'Statistic':>10}  {'Prob':>10}",
            "  " + "-" * 66,
        ]
        labels = {
            "gehan_perm": "Gehan (permutation variance)",
            "gehan_hyper": "Gehan (hypergeometric variance)",
            "logrank": "Logrank",
            "peto_peto": "Peto & Peto",
            "peto_prentice": "Peto & Prentice",
        }
        for key, ts in self.tests.items():
            label = labels.get(key, key)
            lines.append(f"  {label:<40}  {ts.statistic:>10.4f}  {ts.probability:>10.4f}")
        return "\n".join(lines)

    def plot(self, ax: Any = None, labels: tuple[str, str] | None = None, **kwargs: Any) -> Any:
        """Overlay KM curves for both groups.

        Requires matplotlib (``pip install asurv[plot]``).
        """
        from asurv._plots import plot_km

        lbl1, lbl2 = (labels or ("Group I", "Group II"))
        if self.km_group1 is not None:
            ax = plot_km(
                self.km_group1.table,
                self.km_group1.censored_points,
                label=lbl1,
                ax=ax,
                **kwargs,
            )
        if self.km_group2 is not None:
            ax = plot_km(
                self.km_group2.table,
                self.km_group2.censored_points,
                label=lbl2,
                ax=ax,
                **kwargs,
            )
        if ax is not None:
            try:
                ax.legend()
            except Exception:
                pass
        return ax

    def __repr__(self) -> str:
        keys = list(self.tests.keys())
        return (
            f"TwoSampleResult(n={self.n}, n_group1={self.n_group1}, "
            f"n_group2={self.n_group2}, tests={keys})"
        )


def twosample(
    df: pd.DataFrame,
    value: str,
    group: str,
    upper_limit: str,
    groups: tuple[Any, Any],
    title: str = "Two-Sample Test",
    variable_name: str | None = None,
    details: bool = False,
    full_km: bool = True,
    differential: tuple[float, int, float] | None = None,
    print_data: bool = False,
) -> TwoSampleResult:
    """Run the ASURV two-sample tests on a censored dataset.

    Parameters
    ----------
    df:
        DataFrame with value, group, and censoring columns.
    value:
        Column name of the measured values.
    group:
        Column name of the group identifier.
    upper_limit:
        Boolean column: ``True`` = upper limit.
    groups:
        Tuple of exactly two group labels to compare, e.g. ``("Normal", "Starburst")``.
    title:
        Analysis title (≤80 characters).
    variable_name:
        Variable label (≤9 characters; defaults to the value column name).
    details:
        Print full test details in the Fortran output.
    full_km:
        Print the full KM table for each group.
    differential:
        Differential KM histogram parameters ``(start, n_bins, bin_size)``.
    print_data:
        Include input data in the output.

    Returns
    -------
    TwoSampleResult
        Contains five test statistics plus per-group KM results.
    """
    from asurv._backend import run
    from asurv._commands import twost_command_lines, twost_data_text
    from asurv._ind import twost_ind
    from asurv._parsers import parse_twosample
    from asurv._validate import validate_title, validate_twosample_df

    validate_title(title)
    validate_twosample_df(df, value, group, upper_limit, groups)

    vname = (variable_name or value)[:9]
    g1, g2 = groups
    all_groups = sorted(df[group].unique().tolist())

    # Map group labels to contiguous integers starting at 0
    group_map: dict[Any, int] = {g: i for i, g in enumerate(all_groups)}
    g1_id = group_map[g1]
    g2_id = group_map[g2]
    group_ids = list(range(len(all_groups)))

    g1_name = str(g1)[:9]
    g2_name = str(g2)[:9]

    ind = twost_ind(df, upper_limit)
    data_text = twost_data_text(df, value, group, ind, group_map)
    cmd = twost_command_lines(
        data_filename="data.dat",
        title=title,
        nvar=1,
        column=1,
        variable_name=vname,
        group_ids=group_ids,
        first_group_id=g1_id,
        second_group_id=g2_id,
        first_group_name=g1_name,
        second_group_name=g2_name,
        output_filename="twost.out",
        details=details,
        full_km=full_km,
        differential=differential,
        print_data=print_data,
    )
    raw = run("TWOST", data_text, cmd,
              command_filename="twost.com", output_filename="twost.out")
    return parse_twosample(raw)
