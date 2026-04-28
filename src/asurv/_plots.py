"""Matplotlib helpers shared by KMResult, TwoSampleResult, and BivariateResult.

matplotlib is an optional dependency.  All functions in this module raise
ImportError with a helpful message if matplotlib is not installed.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd


def _require_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415
        return plt
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install it with: pip install matplotlib  (or pip install asurv[plot])"
        ) from exc


def plot_km(
    table: "pd.DataFrame",
    censored_points: list[float],
    label: str | None = None,
    ax: Any = None,
    **kwargs: Any,
) -> Any:
    """Draw a KM survival step function with censoring tick marks.

    Parameters
    ----------
    table:
        DataFrame with columns ``from``, ``to``, ``S``, ``error``.
    censored_points:
        X-values of censored observations (drawn as vertical ticks on the curve).
    label:
        Legend label for the step curve.
    ax:
        Existing matplotlib Axes.  If ``None``, uses ``plt.gca()``.
    """
    plt = _require_matplotlib()
    if ax is None:
        ax = plt.gca()

    xs = [row["from"] for _, row in table.iterrows()]
    ys = [row["S"] for _, row in table.iterrows()]
    for _, row in table.iterrows():
        to = row["to"] if row["to"] != float("inf") else row["from"] * 1.1 + 1
        xs.append(to)
        ys.append(row["S"])
    xs_sorted, ys_sorted = zip(*sorted(zip(xs, ys)))
    color = kwargs.pop("color", None)
    line, = ax.step(xs_sorted, ys_sorted, where="post", label=label, color=color, **kwargs)

    for x in censored_points:
        y = _km_at(table, x)
        ax.plot(x, y, "+", color=line.get_color(), markersize=8)

    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("Value")
    ax.set_ylabel("Kaplan-Meier S(x)")
    return ax


def _km_at(table: "pd.DataFrame", x: float) -> float:
    """Interpolate KM survival value at x from the step table."""
    s = 1.0
    for _, row in table.iterrows():
        if x >= row["from"]:
            s = row["S"]
        else:
            break
    return s


def plot_scatter_regression(
    x: "Any",
    y: "Any",
    x_upper: "Any | None" = None,
    y_upper: "Any | None" = None,
    intercept: float | None = None,
    slopes: "list[float] | None" = None,
    bj_intercept: float | None = None,
    bj_slope: float | None = None,
    ax: Any = None,
    x_label: str = "X",
    y_label: str = "Y",
) -> Any:
    """Scatter plot with optional regression line(s) and censoring arrows."""
    import numpy as np

    plt = _require_matplotlib()
    if ax is None:
        ax = plt.gca()

    x = list(x)
    y = list(y)
    n = len(x)
    x_upper_arr = [bool(x_upper[i]) if x_upper is not None else False for i in range(n)]
    y_upper_arr = [bool(y_upper[i]) if y_upper is not None else False for i in range(n)]

    det_x = [x[i] for i in range(n) if not x_upper_arr[i] and not y_upper_arr[i]]
    det_y = [y[i] for i in range(n) if not x_upper_arr[i] and not y_upper_arr[i]]
    ax.scatter(det_x, det_y, zorder=3, label="Detected")

    for i in range(n):
        if y_upper_arr[i]:
            ax.annotate("", xy=(x[i], y[i]), xytext=(x[i], y[i] + 0.3),
                        arrowprops=dict(arrowstyle="->", color="gray"))
        if x_upper_arr[i]:
            ax.annotate("", xy=(x[i], y[i]), xytext=(x[i] + 0.3, y[i]),
                        arrowprops=dict(arrowstyle="->", color="gray", linestyle="--"))

    x_arr = np.array(x)
    x_fit = np.linspace(x_arr.min(), x_arr.max(), 200)

    if intercept is not None and slopes:
        y_fit = intercept + slopes[0] * x_fit
        ax.plot(x_fit, y_fit, "-", label="EM fit")

    if bj_intercept is not None and bj_slope is not None:
        y_fit_bj = bj_intercept + bj_slope * x_fit
        ax.plot(x_fit, y_fit_bj, "--", label="B-J fit")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    return ax
