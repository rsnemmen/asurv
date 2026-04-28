"""End-to-end tests for asurv.km() using the gal1 fixture."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import asurv


def _load_gal1(path: Path) -> pd.DataFrame:
    df = pd.read_fwf(path, widths=[4, 10], header=None, names=["ind", "value"])
    df["upper_limit"] = df["ind"] < 0
    return df


@pytest.mark.integration
def test_km_mean(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies", variable_name="IR",
                 full_km=True, differential=(25.0, 5, 2.0))
    assert abs(r.mean - 27.767) < 0.01


@pytest.mark.integration
def test_km_n(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies")
    assert r.n == 6
    assert r.n_upper_limits == 3


@pytest.mark.integration
def test_km_table(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies")
    assert list(r.table.columns) == ["from", "to", "S", "error"]
    assert len(r.table) >= 3


@pytest.mark.integration
def test_km_censored_points(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies")
    assert len(r.censored_points) == 3


@pytest.mark.integration
def test_km_differential(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies",
                 differential=(25.0, 5, 2.0))
    assert r.differential is not None
    assert len(r.differential) == 5


@pytest.mark.integration
def test_km_result_repr(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies")
    assert "KMResult" in repr(r)


@pytest.mark.integration
def test_km_summary(gal1_dat: Path) -> None:
    df = _load_gal1(gal1_dat)
    r = asurv.km(df, value="value", upper_limit="upper_limit",
                 title="IR Luminosities of Galaxies")
    s = r.summary()
    assert "Mean" in s
    assert "N = 6" in s
