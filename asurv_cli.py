#!/usr/bin/env python3
"""Python CLI shell for driving the legacy ASURV executable."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


LEGACY_NAME_LIMIT = 9
MAX_ROWS = 500
MAX_VARIABLES = 4
MAX_BIVAR_INDEPENDENT = 2
DEFAULT_EXECUTABLE = Path(__file__).resolve().with_name("asurv")

METHOD_TO_NUMBER = {
    "cox": 1,
    "kendall": 2,
    "spearman": 3,
    "em": 4,
    "buckley_james": 5,
    "schmitt": 6,
}

METHOD_ALIASES = {
    "bj": "buckley_james",
    "buckley-james": "buckley_james",
    "buckley_james": "buckley_james",
}


class CliError(Exception):
    """Base exception for CLI failures."""


class ValidationError(CliError):
    """Raised for request validation errors."""


class ExecutionError(CliError):
    """Raised for executable/runtime errors."""


def normalize_method(name: str) -> str:
    normalized = name.strip().lower().replace("-", "_")
    normalized = METHOD_ALIASES.get(normalized, normalized)
    if normalized not in METHOD_TO_NUMBER:
        raise ValidationError(f"unsupported method: {name}")
    return normalized


def validate_title(title: str) -> None:
    if not title.strip():
        raise ValidationError("title must not be blank")
    if len(title) > 80:
        raise ValidationError("title must be 80 characters or fewer")


def validate_legacy_name(value: str, label: str) -> None:
    if not value.strip():
        raise ValidationError(f"{label} must not be blank")
    if len(value) > LEGACY_NAME_LIMIT:
        raise ValidationError(
            f"{label} must be {LEGACY_NAME_LIMIT} characters or fewer: {value!r}"
        )


def validate_positive_int(value: int, label: str) -> None:
    if value <= 0:
        raise ValidationError(f"{label} must be positive")


def require_diff_args(args: argparse.Namespace) -> tuple[int, float, float, int]:
    has_any = any(
        value is not None
        for value in (args.diff_start, args.diff_bins, args.diff_size)
    )
    if not has_any:
        return 0, 0.0, 1, 1.0
    if None in (args.diff_start, args.diff_bins, args.diff_size):
        raise ValidationError(
            "differential KM requires --diff-start, --diff-bins, and --diff-size"
        )
    validate_positive_int(args.diff_bins, "differential KM bin count")
    if args.diff_size <= 0.0:
        raise ValidationError("differential KM bin size must be positive")
    return 1, float(args.diff_start), int(args.diff_bins), float(args.diff_size)


def read_nonblank_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped:
            rows.append(stripped.split())
    return rows


def validate_km_data(path: Path, nvar: int) -> None:
    rows = read_nonblank_rows(path)
    validate_positive_int(len(rows), "number of data rows")
    if len(rows) > MAX_ROWS:
        raise ValidationError(f"data file has more than {MAX_ROWS} rows")
    expected = 2 * nvar
    for index, row in enumerate(rows, start=1):
        if len(row) != expected:
            raise ValidationError(
                f"row {index} in {path} has {len(row)} fields, expected {expected}"
            )


def validate_twost_data(path: Path, nvar: int) -> set[int]:
    rows = read_nonblank_rows(path)
    validate_positive_int(len(rows), "number of data rows")
    if len(rows) > MAX_ROWS:
        raise ValidationError(f"data file has more than {MAX_ROWS} rows")
    expected = 1 + (2 * nvar)
    groups: set[int] = set()
    for index, row in enumerate(rows, start=1):
        if len(row) != expected:
            raise ValidationError(
                f"row {index} in {path} has {len(row)} fields, expected {expected}"
            )
        try:
            groups.add(int(row[0]))
        except ValueError as exc:
            raise ValidationError(
                f"row {index} in {path} has a non-integer group indicator"
            ) from exc
    return groups


def validate_bivar_data(path: Path, n_independent: int) -> tuple[set[int], bool]:
    rows = read_nonblank_rows(path)
    validate_positive_int(len(rows), "number of data rows")
    if len(rows) > MAX_ROWS:
        raise ValidationError(f"data file has more than {MAX_ROWS} rows")
    expected = 2 + n_independent
    indicators: set[int] = set()
    x_or_dual_censoring = False
    for index, row in enumerate(rows, start=1):
        if len(row) != expected:
            raise ValidationError(
                f"row {index} in {path} has {len(row)} fields, expected {expected}"
            )
        try:
            indicator = int(row[0])
        except ValueError as exc:
            raise ValidationError(
                f"row {index} in {path} has a non-integer censoring indicator"
            ) from exc
        indicators.add(indicator)
        if abs(indicator) in {2, 3, 4}:
            x_or_dual_censoring = True
    return indicators, x_or_dual_censoring


def stage_data_file(source: Path, dest: Path) -> None:
    if not source.is_file():
        raise ValidationError(f"data file does not exist: {source}")
    shutil.copyfile(source, dest)


def stage_name(source: Path, fallback: str) -> str:
    name = source.name
    if len(name) <= LEGACY_NAME_LIMIT:
        return name
    return fallback


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def render_double_a9(first: str, second: str) -> str:
    return f"{first:<9}{second if len(second) == 9 else second.rjust(9)}"


def render_km_command(args: argparse.Namespace, data_name: str, output_name: str) -> None:
    kdiff, diff_start, diff_bins, diff_size = require_diff_args(args)
    lines = [
        data_name,
        args.title,
        str(args.nvar),
        str(args.column),
        args.variable_name,
        "1" if args.full_km else "0",
        str(kdiff),
        f"{diff_start:.3f}",
        str(diff_bins),
        f"{diff_size:.3f}",
        "1" if args.print_data else "0",
        output_name,
    ]
    write_lines(args.command_path, lines)


def render_twost_command(
    args: argparse.Namespace, data_name: str, output_name: str
) -> None:
    kdiff, diff_start, diff_bins, diff_size = require_diff_args(args)
    lines = [
        data_name,
        args.title,
        str(args.nvar),
        str(args.column),
        args.variable_name,
        str(len(args.group_ids)),
        "   ".join(str(group) for group in args.group_ids),
        "   ".join(
            [
                str(args.first_group),
                str(args.second_group),
                "1" if args.details else "0",
                "1" if args.full_km else "0",
            ]
        ),
        str(kdiff),
    ]
    if kdiff:
        lines.extend(
            [
                f"{diff_start:.3f}",
                str(diff_bins),
                f"{diff_size:.3f}",
            ]
        )
    lines.extend(
        [
            "1" if args.print_data else "0",
            args.first_group_name,
            args.second_group_name,
            output_name,
        ]
    )
    write_lines(args.command_path, lines)


def render_bivar_command(
    args: argparse.Namespace, data_name: str, output_name: str
) -> None:
    lines = [
        args.title,
        data_name,
        "   ".join(
            [
                str(args.n_independent),
                str(args.column),
                str(len(args.methods)),
            ]
        ),
        "   ".join(str(METHOD_TO_NUMBER[method]) for method in args.methods),
        render_double_a9(args.x_name, args.y_name),
        "1" if args.print_data else "0",
        output_name,
    ]

    for method in args.methods:
        if method == "spearman":
            lines.append("1" if args.spearman_print_ranks else "0")
        elif method == "em":
            coeffs = args.em_initial_coefficients or []
            expected_coeffs = args.n_independent + 2
            if coeffs and len(coeffs) != expected_coeffs:
                raise ValidationError(
                    f"EM initial coefficients require {expected_coeffs} values"
                )
            if not coeffs:
                coeffs = [0.0] * expected_coeffs
            lines.append(f"{args.em_tolerance:.3E}")
            lines.append("".join(f"{value:10.3f}" for value in coeffs))
            lines.append(str(args.em_max_iterations))
        elif method == "buckley_james":
            lines.append(f"{args.bj_tolerance:.3E}")
            lines.append(str(args.bj_max_iterations))
        elif method == "schmitt":
            if args.schmitt_x_bins is None or args.schmitt_y_bins is None:
                raise ValidationError(
                    "Schmitt method requires --schmitt-x-bins and --schmitt-y-bins"
                )
            lines.append(f"{args.schmitt_x_bins:4d}{args.schmitt_y_bins:4d}")
            lines.append("1" if args.schmitt_use_bin_geometry else "0")
            lines.append(f"{args.schmitt_tolerance:.3E}")
            lines.append(str(args.schmitt_max_iterations))
            if args.schmitt_use_bin_geometry:
                required = (
                    args.schmitt_x_bin_size,
                    args.schmitt_y_bin_size,
                    args.schmitt_x_origin,
                    args.schmitt_y_origin,
                )
                if None in required:
                    raise ValidationError(
                        "Schmitt bin geometry requires sizes and origins"
                    )
                if args.schmitt_x_bin_size <= 0.0 or args.schmitt_y_bin_size <= 0.0:
                    raise ValidationError("Schmitt bin sizes must be positive")
                lines.append(
                    f"{args.schmitt_x_bin_size:10.3f}{args.schmitt_y_bin_size:10.3f}"
                )
                lines.append(
                    f"{args.schmitt_x_origin:10.3f}{args.schmitt_y_origin:10.3f}"
                )
            else:
                lines.append(f"{0.0:10.3f}{0.0:10.3f}")
                lines.append(f"{0.0:10.3f}{0.0:10.3f}")
            lines.append("1" if args.schmitt_print_2d_km else "0")
            bootstrap_iterations = args.schmitt_bootstrap_iterations
            if bootstrap_iterations < 0:
                raise ValidationError("Schmitt bootstrap iterations must be nonnegative")
            lines.append(str(bootstrap_iterations))
            if bootstrap_iterations > 0:
                if args.schmitt_random_seed is None:
                    raise ValidationError(
                        "Schmitt bootstrap requires --schmitt-random-seed"
                    )
                if args.schmitt_random_seed >= 0:
                    raise ValidationError(
                        "Schmitt random seed must be negative when bootstrap is enabled"
                    )
                lines.append(str(args.schmitt_random_seed))

    write_lines(args.command_path, lines)


def ensure_bivar_method_compatibility(
    n_independent: int, methods: list[str], has_x_or_dual_censoring: bool
) -> None:
    method_set = set(methods)
    if n_independent > 1:
        allowed = {"cox", "em", "buckley_james"}
        invalid = method_set - allowed
        if invalid:
            raise ValidationError(
                "multivariate BIVAR requests only support cox, em, and buckley-james"
            )
        return
    if has_x_or_dual_censoring:
        allowed = {"kendall", "spearman", "schmitt"}
        invalid = method_set - allowed
        if invalid:
            raise ValidationError(
                "single-variable BIVAR data with X or dual censoring only supports "
                "kendall, spearman, and schmitt"
            )


def run_legacy(
    executable: Path, workspace: Path, menu_script: str, output_name: str
) -> tuple[str, str]:
    executable = executable.resolve()
    if not executable.is_file():
        raise ExecutionError(f"ASURV executable not found: {executable}")
    completed = subprocess.run(
        [str(executable)],
        cwd=workspace,
        input=menu_script,
        text=True,
        capture_output=True,
        check=False,
    )
    output_path = workspace / output_name
    if completed.returncode != 0:
        raise ExecutionError(
            "ASURV exited with a nonzero status:\n"
            f"{completed.stdout}{completed.stderr}"
        )
    if not output_path.is_file():
        raise ExecutionError(
            f"ASURV did not produce the expected output file: {output_path.name}\n"
            f"{completed.stdout}{completed.stderr}"
        )
    return completed.stdout, output_path.read_text()


def build_common_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Python CLI shell for the legacy ASURV executable."
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help="path to the ASURV executable (default: ./asurv)",
    )
    parser.add_argument(
        "--keep-workspace",
        action="store_true",
        help="keep the temporary run directory instead of deleting it",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    km = subparsers.add_parser("km", help="run the Kaplan-Meier workflow")
    add_km_args(km)

    twost = subparsers.add_parser("twost", help="run the two-sample workflow")
    add_twost_args(twost)

    bivar = subparsers.add_parser("bivar", help="run the bivariate workflow")
    add_bivar_args(bivar)
    return parser


def add_km_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--data-file", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--nvar", type=int, required=True)
    parser.add_argument("--column", type=int, required=True)
    parser.add_argument("--variable-name", required=True)
    parser.add_argument("--full-km", action="store_true")
    parser.add_argument("--diff-start", type=float)
    parser.add_argument("--diff-bins", type=int)
    parser.add_argument("--diff-size", type=float)
    parser.add_argument("--print-data", action="store_true")
    parser.add_argument("--report", type=Path, help="write report text to this path")


def add_twost_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--data-file", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--nvar", type=int, required=True)
    parser.add_argument("--column", type=int, required=True)
    parser.add_argument("--variable-name", required=True)
    parser.add_argument("--group-id", type=int, action="append", dest="group_ids")
    parser.add_argument("--first-group", type=int, required=True)
    parser.add_argument("--second-group", type=int, required=True)
    parser.add_argument("--first-group-name", required=True)
    parser.add_argument("--second-group-name", required=True)
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--full-km", action="store_true")
    parser.add_argument("--diff-start", type=float)
    parser.add_argument("--diff-bins", type=int)
    parser.add_argument("--diff-size", type=float)
    parser.add_argument("--print-data", action="store_true")
    parser.add_argument("--report", type=Path, help="write report text to this path")


def add_bivar_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--data-file", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--n-independent", type=int, required=True)
    parser.add_argument("--column", type=int, default=1)
    parser.add_argument("--x-name", required=True)
    parser.add_argument("--y-name", required=True)
    parser.add_argument(
        "--method",
        action="append",
        required=True,
        help="repeat for each method: cox, kendall, spearman, em, buckley-james, schmitt",
    )
    parser.add_argument("--print-data", action="store_true")
    parser.add_argument("--report", type=Path, help="write report text to this path")
    parser.add_argument("--spearman-print-ranks", action="store_true")
    parser.add_argument("--em-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--em-max-iterations", type=int, default=50)
    parser.add_argument("--em-initial-coefficients", type=float, nargs="+")
    parser.add_argument("--bj-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--bj-max-iterations", type=int, default=50)
    parser.add_argument("--schmitt-x-bins", type=int)
    parser.add_argument("--schmitt-y-bins", type=int)
    parser.add_argument("--schmitt-use-bin-geometry", action="store_true")
    parser.add_argument("--schmitt-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--schmitt-max-iterations", type=int, default=50)
    parser.add_argument("--schmitt-x-bin-size", type=float)
    parser.add_argument("--schmitt-y-bin-size", type=float)
    parser.add_argument("--schmitt-x-origin", type=float)
    parser.add_argument("--schmitt-y-origin", type=float)
    parser.add_argument("--schmitt-print-2d-km", action="store_true")
    parser.add_argument("--schmitt-bootstrap-iterations", type=int, default=0)
    parser.add_argument("--schmitt-random-seed", type=int)


def validate_common_dimensions(nvar: int, column: int) -> None:
    validate_positive_int(nvar, "number of variables")
    if nvar > MAX_VARIABLES:
        raise ValidationError(
            f"number of variables must not exceed the current ASURV limit ({MAX_VARIABLES})"
        )
    if column < 1 or column > nvar:
        raise ValidationError("selected variable index is out of range")


def validate_bivar_dimensions(n_independent: int, column: int) -> None:
    validate_positive_int(n_independent, "number of independent variables")
    if n_independent > MAX_BIVAR_INDEPENDENT:
        raise ValidationError(
            "number of independent variables must not exceed the current "
            f"ASURV bivariate limit ({MAX_BIVAR_INDEPENDENT})"
        )
    if column < 1 or column > n_independent:
        raise ValidationError("selected variable index is out of range")


def run_km(args: argparse.Namespace) -> int:
    validate_title(args.title)
    validate_common_dimensions(args.nvar, args.column)
    validate_legacy_name(args.variable_name, "variable name")
    validate_km_data(args.data_file, args.nvar)
    temp_dir = Path(tempfile.mkdtemp())
    try:
        workspace = temp_dir
        data_name = stage_name(args.data_file, "km.dat")
        stage_data_file(args.data_file, workspace / data_name)
        args.command_path = workspace / "km.com"
        render_km_command(args, data_name, "km.out")
        stdout_text, report_text = run_legacy(
            args.executable,
            workspace,
            "\n\n1\n1\ny\nkm.com\nn\n",
            "km.out",
        )
        emit_results(args, workspace, report_text, stdout_text)
        if args.keep_workspace:
            print(f"workspace kept at {workspace}", file=sys.stderr)
        else:
            shutil.rmtree(workspace)
    except Exception:
        if not args.keep_workspace:
            shutil.rmtree(temp_dir, ignore_errors=True)
        raise
    return 0


def run_twost(args: argparse.Namespace) -> int:
    validate_title(args.title)
    validate_common_dimensions(args.nvar, args.column)
    validate_legacy_name(args.variable_name, "variable name")
    validate_legacy_name(args.first_group_name, "first group name")
    validate_legacy_name(args.second_group_name, "second group name")
    data_groups = validate_twost_data(args.data_file, args.nvar)
    if not args.group_ids:
        raise ValidationError("at least one --group-id is required")
    if len(set(args.group_ids)) != len(args.group_ids):
        raise ValidationError("group ids must be unique")
    if args.first_group == args.second_group:
        raise ValidationError("first and second group must differ")
    if args.first_group not in args.group_ids or args.second_group not in args.group_ids:
        raise ValidationError("selected groups must appear in --group-id")
    if not data_groups.issubset(set(args.group_ids)):
        raise ValidationError("data file contains group ids not listed with --group-id")
    temp_dir = Path(tempfile.mkdtemp())
    try:
        workspace = temp_dir
        data_name = stage_name(args.data_file, "ts.dat")
        stage_data_file(args.data_file, workspace / data_name)
        args.command_path = workspace / "ts.com"
        render_twost_command(args, data_name, "ts.out")
        stdout_text, report_text = run_legacy(
            args.executable,
            workspace,
            "\n\n1\n2\ny\nts.com\nn\n",
            "ts.out",
        )
        emit_results(args, workspace, report_text, stdout_text)
        if args.keep_workspace:
            print(f"workspace kept at {workspace}", file=sys.stderr)
        else:
            shutil.rmtree(workspace)
    except Exception:
        if not args.keep_workspace:
            shutil.rmtree(temp_dir, ignore_errors=True)
        raise
    return 0


def run_bivar(args: argparse.Namespace) -> int:
    validate_title(args.title)
    validate_bivar_dimensions(args.n_independent, args.column)
    validate_legacy_name(args.x_name, "independent variable name")
    validate_legacy_name(args.y_name, "dependent variable name")
    args.methods = [normalize_method(method) for method in args.method]
    if len(set(args.methods)) != len(args.methods):
        raise ValidationError("methods must not be repeated")
    _, has_x_or_dual_censoring = validate_bivar_data(args.data_file, args.n_independent)
    ensure_bivar_method_compatibility(
        args.n_independent,
        args.methods,
        has_x_or_dual_censoring,
    )
    if args.em_tolerance <= 0.0 or args.bj_tolerance <= 0.0 or args.schmitt_tolerance <= 0.0:
        raise ValidationError("tolerance values must be positive")
    validate_positive_int(args.em_max_iterations, "EM max iterations")
    validate_positive_int(args.bj_max_iterations, "Buckley-James max iterations")
    validate_positive_int(args.schmitt_max_iterations, "Schmitt max iterations")
    temp_dir = Path(tempfile.mkdtemp())
    try:
        workspace = temp_dir
        data_name = stage_name(args.data_file, "bv.dat")
        stage_data_file(args.data_file, workspace / data_name)
        args.command_path = workspace / "bv.com"
        render_bivar_command(args, data_name, "bv.out")
        stdout_text, report_text = run_legacy(
            args.executable,
            workspace,
            "\n\n2\ny\nbv.com\nn\n",
            "bv.out",
        )
        emit_results(args, workspace, report_text, stdout_text)
        if args.keep_workspace:
            print(f"workspace kept at {workspace}", file=sys.stderr)
        else:
            shutil.rmtree(workspace)
    except Exception:
        if not args.keep_workspace:
            shutil.rmtree(temp_dir, ignore_errors=True)
        raise
    return 0


def emit_results(
    args: argparse.Namespace, workspace: Path, report_text: str, stdout_text: str
) -> None:
    if args.report is not None:
        args.report.write_text(report_text)
    else:
        sys.stdout.write(report_text)
    if args.keep_workspace:
        debug_log = workspace / "stdout.log"
        debug_log.write_text(stdout_text)


def main(argv: list[str] | None = None) -> int:
    parser = build_common_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "km":
            return run_km(args)
        if args.command == "twost":
            return run_twost(args)
        if args.command == "bivar":
            return run_bivar(args)
    except CliError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
