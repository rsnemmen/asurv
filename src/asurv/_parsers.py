"""Parse ASURV .out files into structured result objects.

All parsers accept the raw output text and return the corresponding result
dataclass.  They raise AsurvParseError (with the raw text attached) when
an expected section is missing.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from asurv.km import KMResult
    from asurv.twosample import TestStat, TwoSampleResult
    from asurv.bivar import (
        BHKResult,
        BJResult,
        BivariateResult,
        CoxResult,
        EMResult,
        SpearmanResult,
        TwoKMResult,
    )


def _require(match: re.Match | None, pattern: str, raw: str) -> re.Match:
    from asurv import AsurvParseError

    if match is None:
        raise AsurvParseError(
            f"could not find expected pattern in ASURV output: {pattern!r}",
            raw_output=raw,
        )
    return match


# ---------------------------------------------------------------------------
# KM parser
# ---------------------------------------------------------------------------

def parse_km(text: str) -> "KMResult":
    from asurv import AsurvParseError
    from asurv.km import KMResult

    raw = text

    # header counts
    m = _require(
        re.search(r"# OF DATA POINTS\s*:\s*(\d+)\s*# OF UPPER LIMITS\s*:\s*(\d+)", raw),
        "# OF DATA POINTS",
        raw,
    )
    n = int(m.group(1))
    n_upper = int(m.group(2))

    # KM table
    table_rows: list[dict] = []
    censored_points: list[float] = []
    for line in raw.splitlines():
        # detected step: "FROM   26.900   TO   28.500       0.375       0.213"
        mrow = re.match(
            r"\s*FROM\s+([\d.]+)\s+TO\s+([\d.]+)\s+([\d.]+)\s*([\d.]+)?", line
        )
        if mrow:
            table_rows.append({
                "from": float(mrow.group(1)),
                "to": float(mrow.group(2)),
                "S": float(mrow.group(3)),
                "error": float(mrow.group(4)) if mrow.group(4) else float("nan"),
            })
            continue
        # final step: "FROM   30.100   ONWARDS           0.000       0.000"
        monwards = re.match(
            r"\s*FROM\s+([\d.]+)\s+ONWARDS\s+([\d.]+)\s*([\d.]+)?", line
        )
        if monwards:
            table_rows.append({
                "from": float(monwards.group(1)),
                "to": float("inf"),
                "S": float(monwards.group(2)),
                "error": float(monwards.group(3)) if monwards.group(3) else float("nan"),
            })
            continue
        # censored point: "       27.600 C "
        mc = re.match(r"\s*([\d.]+)\s+C\s*$", line)
        if mc:
            censored_points.append(float(mc.group(1)))

    if not table_rows:
        raise AsurvParseError("KM table not found in output", raw_output=raw)

    table = pd.DataFrame(table_rows)

    # percentiles
    percentiles: dict[int, float] | None = None
    mp = re.search(
        r"PERCENTILES\s+(\d+)\s+TH\s+(\d+)\s+TH\s+(\d+)\s+TH\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        raw,
        re.DOTALL,
    )
    if mp:
        percentiles = {
            int(mp.group(1)): float(mp.group(4)),
            int(mp.group(2)): float(mp.group(5)),
            int(mp.group(3)): float(mp.group(6)),
        }

    # mean
    mm = _require(re.search(r"MEAN=\s*([\d.]+)\s*\+/-\s*([\d.]+)", raw), "MEAN=", raw)
    mean = float(mm.group(1))
    mean_err = float(mm.group(2))

    # differential KM
    differential: pd.DataFrame | None = None
    if "DIFFERENTIAL KM ESTIMATOR" in raw:
        diff_rows: list[dict] = []
        in_diff = False
        for line in raw.splitlines():
            if "DIFFERENTIAL KM ESTIMATOR" in line:
                in_diff = True
                continue
            if in_diff:
                md = re.match(r"\s*([\d.]+)\s+([\d.]+)\s*$", line)
                if md:
                    diff_rows.append({
                        "bin_center": float(md.group(1)),
                        "D": float(md.group(2)),
                    })
                elif re.match(r"\s*-+", line):
                    break
        if diff_rows:
            differential = pd.DataFrame(diff_rows)

    return KMResult(
        n=n,
        n_upper_limits=n_upper,
        table=table,
        censored_points=censored_points,
        percentiles=percentiles,
        mean=mean,
        mean_err=mean_err,
        differential=differential,
        raw_output=raw,
    )


# ---------------------------------------------------------------------------
# Two-sample parser
# ---------------------------------------------------------------------------

_TEST_HEADERS = {
    "gehan_perm": r"GEHAN`S GENERALIZED WILCOXON TEST -- PERMUTATION VARIANCE",
    "gehan_hyper": r"GEHAN`S GENERALIZED WILCOXON TEST -- HYPERGEOMETRIC VARIANCE",
    "logrank": r"LOGRANK TEST",
    "peto_peto": r"PETO & PETO GENERALIZED WILCOXON TEST",
    "peto_prentice": r"PETO & PRENTICE GENERALIZED WILCOXON TEST",
}


def parse_twosample(text: str) -> "TwoSampleResult":
    from asurv import AsurvParseError
    from asurv.twosample import TestStat, TwoSampleResult

    raw = text

    # header counts
    m = _require(
        re.search(
            r"# OF DATA POINTS\s*:\s*(\d+),\s*# OF UPPER LIMITS\s*:\s*(\d+)", raw
        ),
        "# OF DATA POINTS",
        raw,
    )
    n = int(m.group(1))
    n_upper = int(m.group(2))

    mg1 = _require(
        re.search(r"# OF DATA POINTS IN GROUP I\s*:\s*(\d+)", raw),
        "GROUP I count",
        raw,
    )
    mg2 = _require(
        re.search(r"# OF DATA POINTS IN GROUP II\s*:\s*(\d+)", raw),
        "GROUP II count",
        raw,
    )
    n1 = int(mg1.group(1))
    n2 = int(mg2.group(1))

    # five test statistics
    tests: dict[str, "TestStat"] = {}
    for key, header in _TEST_HEADERS.items():
        ms = re.search(
            rf"{re.escape(header)}.*?TEST STATISTIC\s*=\s*([\d.]+).*?PROBABILITY\s*=\s*([\d.]+)",
            raw,
            re.DOTALL,
        )
        if ms:
            tests[key] = TestStat(
                statistic=float(ms.group(1)),
                probability=float(ms.group(2)),
            )

    if not tests:
        raise AsurvParseError("no test statistics found in two-sample output", raw_output=raw)

    # per-group KM tables — they follow the test blocks
    km_sections = list(re.finditer(r"KAPLAN-MEIER ESTIMATOR", raw))
    km1: "KMResult | None" = None
    km2: "KMResult | None" = None
    if len(km_sections) >= 1:
        start1 = km_sections[0].start()
        end1 = km_sections[1].start() if len(km_sections) >= 2 else len(raw)
        try:
            km1 = parse_km(raw[start1:end1])
        except AsurvParseError:
            km1 = None
    if len(km_sections) >= 2:
        start2 = km_sections[1].start()
        try:
            km2 = parse_km(raw[start2:])
        except AsurvParseError:
            km2 = None

    return TwoSampleResult(
        n=n,
        n_upper_limits=n_upper,
        n_group1=n1,
        n_group2=n2,
        tests=tests,
        km_group1=km1,
        km_group2=km2,
        raw_output=raw,
    )


# ---------------------------------------------------------------------------
# Bivariate parser
# ---------------------------------------------------------------------------

def _parse_cox(text: str) -> "CoxResult | None":
    from asurv.bivar import CoxResult

    m = re.search(
        r"CORRELATION TEST BY COX PROPORTIONAL HAZARD MODEL.*?"
        r"GLOBAL CHI SQUARE\s*=\s*([\d.]+)\s*WITH\s*(\d+)\s*DEGREES.*?"
        r"PROBABILITY\s*=\s*([\d.]+)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None
    return CoxResult(
        chi2=float(m.group(1)),
        dof=int(m.group(2)),
        probability=float(m.group(3)),
    )


def _parse_bhk(text: str) -> "BHKResult | None":
    from asurv.bivar import BHKResult

    m = re.search(
        r"GENERALIZED KENDALL.*?Z\s*=\s*([-\d.]+).*?PROBABILITY\s*=\s*([\d.]+)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None
    return BHKResult(z=float(m.group(1)), probability=float(m.group(2)))


def _parse_spearman(text: str) -> "SpearmanResult | None":
    from asurv.bivar import SpearmanResult

    m = re.search(
        r"GENERALIZED SPEARMAN.*?RHO\s*=\s*([-\d.]+).*?Z\s*=\s*([-\d.]+).*?PROBABILITY\s*=\s*([\d.]+)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None
    return SpearmanResult(
        rho=float(m.group(1)),
        z=float(m.group(2)),
        probability=float(m.group(3)),
    )


def _parse_em(text: str) -> "EMResult | None":
    from asurv.bivar import EMResult

    m = re.search(
        r"LINEAR REGRESSION BY PARAMETRIC EM ALGORITHM.*?"
        r"INTERCEPT COEFF\s*:\s*([-\d.]+)\s*\+/-\s*([-\d.]+)(.*?)"
        r"STANDARD DEVIATION\s*:\s*([\d.]+).*?"
        r"ITERATIONS\s*:\s*(\d+)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None
    intercept = float(m.group(1))
    intercept_err = float(m.group(2))
    slopes_text = m.group(3)
    sigma = float(m.group(4))
    iterations = int(m.group(5))

    slopes: list[tuple[float, float]] = []
    for sm in re.finditer(
        r"SLOPE COEFF\s+\d+\s*:\s*([-\d.]+)\s*\+/-\s*([-\d.]+)", slopes_text
    ):
        slopes.append((float(sm.group(1)), float(sm.group(2))))

    return EMResult(
        intercept=intercept,
        intercept_err=intercept_err,
        slopes=slopes,
        sigma=sigma,
        iterations=iterations,
    )


def _parse_bj(text: str) -> "BJResult | None":
    from asurv.bivar import BJResult

    m = re.search(
        r"BUCKLEY.JAMES.*?"
        r"INTERCEPT\s*:\s*([-\d.]+)\s*\+/-\s*([-\d.]+).*?"
        r"SLOPE\s*:\s*([-\d.]+)\s*\+/-\s*([-\d.]+)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None
    return BJResult(
        intercept=float(m.group(1)),
        intercept_err=float(m.group(2)),
        slope=float(m.group(3)),
        slope_err=float(m.group(4)),
    )


def _parse_twokm(text: str) -> "TwoKMResult | None":
    from asurv.bivar import TwoKMResult

    if "SCHMITT" not in text.upper() and "TWOKM" not in text.upper():
        return None
    m = re.search(
        r"SLOPE\s*=\s*([-\d.]+).*?INTERCEPT\s*=\s*([-\d.]+)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None
    return TwoKMResult(
        slope=float(m.group(1)),
        intercept=float(m.group(2)),
    )


def parse_bivar(text: str) -> "BivariateResult":
    from asurv import AsurvParseError
    from asurv.bivar import BivariateResult

    raw = text

    # header counts
    mn = re.search(r"NUMBER OF DATA POINTS\s*:\s*(\d+)", raw)
    if mn is None:
        raise AsurvParseError("bivariate output: NUMBER OF DATA POINTS not found", raw_output=raw)
    n = int(mn.group(1))

    upper_limits: dict[str, int] = {"y": 0, "x": 0, "both": 0, "mix": 0}
    mul = re.search(
        r"UPPER LIMITS IN\s+Y\s+X\s+BOTH\s+MIX\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)",
        raw,
    )
    if mul:
        upper_limits = {
            "y": int(mul.group(1)),
            "x": int(mul.group(2)),
            "both": int(mul.group(3)),
            "mix": int(mul.group(4)),
        }

    return BivariateResult(
        n=n,
        upper_limits=upper_limits,
        cox=_parse_cox(raw),
        bhk=_parse_bhk(raw),
        spearman=_parse_spearman(raw),
        em=_parse_em(raw),
        bj=_parse_bj(raw),
        twokm=_parse_twokm(raw),
        raw_output=raw,
    )
