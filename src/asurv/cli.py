"""Command-line interface for the ASURV Python package.

This module re-implements the argparse interface from the old asurv_cli.py,
wired to the new asurv Python API.  The legacy asurv_cli.py is kept as a
thin shim that calls main() here for backward compatibility.

Entry point: ``asurv-py`` (see pyproject.toml [project.scripts]).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_fixed_width_dat(path: Path, widths: list[int]) -> pd.DataFrame:
    """Read an ASURV fixed-width data file into a DataFrame."""
    return pd.read_fwf(path, widths=widths, header=None)


# ---------------------------------------------------------------------------
# km subcommand
# ---------------------------------------------------------------------------

def _add_km_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--data-file", type=Path, required=True,
                   help="ASURV-format data file (IND, value per row)")
    p.add_argument("--title", required=True, help="analysis title")
    p.add_argument("--variable-name", default=None,
                   help="variable label (default: 'value')")
    p.add_argument("--full-km", action="store_true",
                   help="print the full KM table")
    p.add_argument("--diff-start", type=float)
    p.add_argument("--diff-bins", type=int)
    p.add_argument("--diff-size", type=float)
    p.add_argument("--print-data", action="store_true")
    p.add_argument("--report", type=Path,
                   help="write raw output text to this file instead of stdout")
    # legacy compat flags
    p.add_argument("--nvar", type=int, default=1, help=argparse.SUPPRESS)
    p.add_argument("--column", type=int, default=1, help=argparse.SUPPRESS)


def _run_km(args: argparse.Namespace) -> int:
    from asurv import km as asurv_km

    df = _read_fixed_width_dat(args.data_file, [4, 10])
    df.columns = ["ind", "value"]
    df["upper_limit"] = df["ind"] < 0

    diff = None
    if args.diff_start is not None:
        diff = (args.diff_start, args.diff_bins, args.diff_size)

    result = asurv_km(
        df,
        value="value",
        upper_limit="upper_limit",
        title=args.title,
        variable_name=args.variable_name or "value",
        full_km=args.full_km,
        differential=diff,
        print_data=args.print_data,
    )
    _emit(args, result)
    return 0


# ---------------------------------------------------------------------------
# twost subcommand
# ---------------------------------------------------------------------------

def _add_twost_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--data-file", type=Path, required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--variable-name", default=None)
    p.add_argument("--first-group-name", required=True)
    p.add_argument("--second-group-name", required=True)
    p.add_argument("--details", action="store_true")
    p.add_argument("--full-km", action="store_true")
    p.add_argument("--diff-start", type=float)
    p.add_argument("--diff-bins", type=int)
    p.add_argument("--diff-size", type=float)
    p.add_argument("--print-data", action="store_true")
    p.add_argument("--report", type=Path)
    # legacy compat
    p.add_argument("--nvar", type=int, default=1, help=argparse.SUPPRESS)
    p.add_argument("--column", type=int, default=1, help=argparse.SUPPRESS)
    p.add_argument("--group-id", type=int, action="append", dest="group_ids",
                   help=argparse.SUPPRESS)
    p.add_argument("--first-group", type=int, help=argparse.SUPPRESS)
    p.add_argument("--second-group", type=int, help=argparse.SUPPRESS)


def _run_twost(args: argparse.Namespace) -> int:
    from asurv import twosample as asurv_twost

    df = _read_fixed_width_dat(args.data_file, [4, 4, 10])
    df.columns = ["group", "ind", "value"]
    df["upper_limit"] = df["ind"] < 0

    g1 = args.first_group_name
    g2 = args.second_group_name
    all_groups = sorted(df["group"].unique().tolist())

    diff = None
    if args.diff_start is not None:
        diff = (args.diff_start, args.diff_bins, args.diff_size)

    # We use the integer group IDs as group labels in the DataFrame
    result = asurv_twost(
        df,
        value="value",
        group="group",
        upper_limit="upper_limit",
        groups=(all_groups[0], all_groups[1]),
        title=args.title,
        variable_name=args.variable_name or "value",
        details=args.details,
        full_km=args.full_km,
        differential=diff,
        print_data=args.print_data,
    )
    _emit(args, result)
    return 0


# ---------------------------------------------------------------------------
# bivar subcommand
# ---------------------------------------------------------------------------

def _add_bivar_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--data-file", type=Path, required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--x-name", required=True)
    p.add_argument("--y-name", required=True)
    p.add_argument("--method", action="append", required=True, dest="methods",
                   help="repeat for each method: cox, kendall, spearman, em, buckley-james, schmitt")
    p.add_argument("--print-data", action="store_true")
    p.add_argument("--report", type=Path)
    p.add_argument("--spearman-print-ranks", action="store_true")
    p.add_argument("--em-tolerance", type=float, default=1e-5)
    p.add_argument("--em-max-iterations", type=int, default=50)
    p.add_argument("--em-initial-coefficients", type=float, nargs="+")
    p.add_argument("--bj-tolerance", type=float, default=1e-5)
    p.add_argument("--bj-max-iterations", type=int, default=50)
    p.add_argument("--schmitt-x-bins", type=int)
    p.add_argument("--schmitt-y-bins", type=int)
    p.add_argument("--schmitt-use-bin-geometry", action="store_true")
    p.add_argument("--schmitt-tolerance", type=float, default=1e-5)
    p.add_argument("--schmitt-max-iterations", type=int, default=50)
    p.add_argument("--schmitt-x-bin-size", type=float)
    p.add_argument("--schmitt-y-bin-size", type=float)
    p.add_argument("--schmitt-x-origin", type=float)
    p.add_argument("--schmitt-y-origin", type=float)
    p.add_argument("--schmitt-print-2d-km", action="store_true")
    p.add_argument("--schmitt-bootstrap-iterations", type=int, default=0)
    p.add_argument("--schmitt-random-seed", type=int)
    # legacy compat
    p.add_argument("--n-independent", type=int, default=1, help=argparse.SUPPRESS)
    p.add_argument("--column", type=int, default=1, help=argparse.SUPPRESS)


def _run_bivar(args: argparse.Namespace) -> int:
    from asurv import bivar as asurv_bivar

    df = _read_fixed_width_dat(args.data_file, [4, 10, 10])
    df.columns = ["ind", "x", "y"]
    # IND codes: -1 = Y upper limit, -2 = X upper limit, -4 = both
    df["y_upper"] = df["ind"].isin([-1, -4])
    df["x_upper"] = df["ind"].isin([-2, -4])

    result = asurv_bivar(
        df,
        x="x",
        y="y",
        methods=args.methods,
        y_upper="y_upper" if df["y_upper"].any() else None,
        x_upper="x_upper" if df["x_upper"].any() else None,
        title=args.title,
        x_name=args.x_name,
        y_name=args.y_name,
        print_data=args.print_data,
        spearman_print_ranks=args.spearman_print_ranks,
        em_tolerance=args.em_tolerance,
        em_max_iterations=args.em_max_iterations,
        em_initial_coefficients=args.em_initial_coefficients,
        bj_tolerance=args.bj_tolerance,
        bj_max_iterations=args.bj_max_iterations,
        schmitt_x_bins=args.schmitt_x_bins,
        schmitt_y_bins=args.schmitt_y_bins,
        schmitt_use_bin_geometry=args.schmitt_use_bin_geometry,
        schmitt_tolerance=args.schmitt_tolerance,
        schmitt_max_iterations=args.schmitt_max_iterations,
        schmitt_x_bin_size=args.schmitt_x_bin_size,
        schmitt_y_bin_size=args.schmitt_y_bin_size,
        schmitt_x_origin=args.schmitt_x_origin,
        schmitt_y_origin=args.schmitt_y_origin,
        schmitt_print_2d_km=args.schmitt_print_2d_km,
        schmitt_bootstrap_iterations=args.schmitt_bootstrap_iterations,
        schmitt_random_seed=args.schmitt_random_seed,
    )
    _emit(args, result)
    return 0


# ---------------------------------------------------------------------------
# Output helper
# ---------------------------------------------------------------------------

def _emit(args: argparse.Namespace, result: object) -> None:
    text = result.raw_output  # type: ignore[attr-defined]
    if args.report is not None:
        args.report.write_text(text)
    else:
        sys.stdout.write(text)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="asurv-py",
        description="Python interface to the ASURV survival analysis package.",
    )
    subs = parser.add_subparsers(dest="command", required=True)

    _add_km_args(subs.add_parser("km", help="Kaplan-Meier estimator"))
    _add_twost_args(subs.add_parser("twost", help="two-sample tests"))
    _add_bivar_args(subs.add_parser("bivar", help="bivariate correlation and regression"))

    args = parser.parse_args(argv)

    from asurv import AsurvError

    try:
        if args.command == "km":
            return _run_km(args)
        if args.command == "twost":
            return _run_twost(args)
        if args.command == "bivar":
            return _run_bivar(args)
    except AsurvError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
