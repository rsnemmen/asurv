"""End-to-end tests for asurv.bivar() using the gal3 fixture."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import asurv


def _load_gal3(path: Path) -> pd.DataFrame:
    df = pd.read_fwf(path, widths=[4, 10, 10], header=None,
                     names=["ind", "x", "y"])
    df["y_upper"] = df["ind"].isin([-1, -4])
    df["x_upper"] = df["ind"].isin([-2, -4])
    return df


@pytest.mark.integration
def test_bivar_n(gal3_dat: Path) -> None:
    df = _load_gal3(gal3_dat)
    r = asurv.bivar(
        df, x="x", y="y",
        methods=["cox", "em"],
        y_upper="y_upper",
        title="Optical-Infrared Relation",
        x_name="Optical",
        y_name="Infrared",
        em_initial_coefficients=[0.0, 0.0, 0.0],
    )
    assert r.n == 16


@pytest.mark.integration
def test_bivar_cox(gal3_dat: Path) -> None:
    df = _load_gal3(gal3_dat)
    r = asurv.bivar(
        df, x="x", y="y",
        methods=["cox", "em"],
        y_upper="y_upper",
        title="Optical-Infrared Relation",
        x_name="Optical",
        y_name="Infrared",
        em_initial_coefficients=[0.0, 0.0, 0.0],
    )
    assert r.cox is not None
    assert abs(r.cox.chi2 - 18.458) < 0.01
    assert r.cox.dof == 1
    assert r.cox.probability < 0.001


@pytest.mark.integration
def test_bivar_em(gal3_dat: Path) -> None:
    df = _load_gal3(gal3_dat)
    r = asurv.bivar(
        df, x="x", y="y",
        methods=["cox", "em"],
        y_upper="y_upper",
        title="Optical-Infrared Relation",
        x_name="Optical",
        y_name="Infrared",
        em_initial_coefficients=[0.0, 0.0, 0.0],
    )
    assert r.em is not None
    assert abs(r.em.intercept - 0.1703) < 0.01
    assert len(r.em.slopes) == 1
    slope, _ = r.em.slopes[0]
    assert abs(slope - 1.0519) < 0.01
    assert abs(r.em.sigma - 0.3687) < 0.01
    assert r.em.iterations == 27


@pytest.mark.integration
def test_bivar_upper_limits(gal3_dat: Path) -> None:
    df = _load_gal3(gal3_dat)
    r = asurv.bivar(
        df, x="x", y="y",
        methods=["cox"],
        y_upper="y_upper",
        title="Optical-Infrared Relation",
        x_name="Optical",
        y_name="Infrared",
    )
    assert r.upper_limits["y"] == 6
    assert r.upper_limits["x"] == 0


@pytest.mark.integration
def test_bivar_summary(gal3_dat: Path) -> None:
    df = _load_gal3(gal3_dat)
    r = asurv.bivar(
        df, x="x", y="y",
        methods=["cox", "em"],
        y_upper="y_upper",
        title="Optical-Infrared Relation",
        x_name="Optical",
        y_name="Infrared",
        em_initial_coefficients=[0.0, 0.0, 0.0],
    )
    s = r.summary()
    assert "Cox" in s
    assert "EM" in s
