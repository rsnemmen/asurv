"""End-to-end tests for asurv.twosample() using the gal2 fixture."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import asurv


def _load_gal2(path: Path) -> pd.DataFrame:
    df = pd.read_fwf(path, widths=[4, 4, 10], header=None,
                     names=["group", "ind", "value"])
    df["upper_limit"] = df["ind"] < 0
    return df


@pytest.mark.integration
def test_twosample_counts(gal2_dat: Path) -> None:
    df = _load_gal2(gal2_dat)
    r = asurv.twosample(
        df, value="value", group="group", upper_limit="upper_limit",
        groups=(0, 1), title="IR Luminosities of Galaxies",
    )
    assert r.n == 12
    assert r.n_upper_limits == 5
    assert r.n_group1 == 6
    assert r.n_group2 == 6


@pytest.mark.integration
def test_twosample_five_tests(gal2_dat: Path) -> None:
    df = _load_gal2(gal2_dat)
    r = asurv.twosample(
        df, value="value", group="group", upper_limit="upper_limit",
        groups=(0, 1), title="IR Luminosities of Galaxies",
    )
    for key in ("gehan_perm", "gehan_hyper", "logrank", "peto_peto", "peto_prentice"):
        assert key in r.tests


@pytest.mark.integration
def test_twosample_gehan_perm(gal2_dat: Path) -> None:
    df = _load_gal2(gal2_dat)
    r = asurv.twosample(
        df, value="value", group="group", upper_limit="upper_limit",
        groups=(0, 1), title="IR Luminosities of Galaxies",
    )
    ts = r.tests["gehan_perm"]
    assert abs(ts.statistic - 1.652) < 0.01
    assert abs(ts.probability - 0.0986) < 0.001


@pytest.mark.integration
def test_twosample_logrank(gal2_dat: Path) -> None:
    df = _load_gal2(gal2_dat)
    r = asurv.twosample(
        df, value="value", group="group", upper_limit="upper_limit",
        groups=(0, 1), title="IR Luminosities of Galaxies",
    )
    ts = r.tests["logrank"]
    assert abs(ts.statistic - 1.814) < 0.01
    assert abs(ts.probability - 0.0696) < 0.001


@pytest.mark.integration
def test_twosample_km_groups(gal2_dat: Path) -> None:
    df = _load_gal2(gal2_dat)
    r = asurv.twosample(
        df, value="value", group="group", upper_limit="upper_limit",
        groups=(0, 1), title="IR Luminosities of Galaxies",
    )
    assert r.km_group1 is not None
    assert r.km_group2 is not None


@pytest.mark.integration
def test_twosample_summary(gal2_dat: Path) -> None:
    df = _load_gal2(gal2_dat)
    r = asurv.twosample(
        df, value="value", group="group", upper_limit="upper_limit",
        groups=(0, 1), title="IR Luminosities of Galaxies",
    )
    s = r.summary()
    assert "Logrank" in s
    assert "Gehan" in s
