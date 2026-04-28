"""Unit tests for _parsers.py using canned output text from asurv.etc.

These tests do not invoke the Fortran binary — they parse fixed strings and
verify the regex logic is correct.
"""

from __future__ import annotations

import math

import pytest

from asurv._parsers import parse_bivar, parse_km, parse_twosample


class TestParseKM:
    def test_basic_fields(self, gal1_out: str) -> None:
        r = parse_km(gal1_out)
        assert r.n == 6
        assert r.n_upper_limits == 3
        assert abs(r.mean - 27.767) < 0.001
        assert abs(r.mean_err - 0.515) < 0.001

    def test_table_shape(self, gal1_out: str) -> None:
        r = parse_km(gal1_out)
        assert list(r.table.columns) == ["from", "to", "S", "error"]
        assert len(r.table) >= 3

    def test_first_row(self, gal1_out: str) -> None:
        r = parse_km(gal1_out)
        row = r.table.iloc[0]
        assert row["from"] == pytest.approx(0.0)
        assert row["to"] == pytest.approx(26.9)
        assert row["S"] == pytest.approx(1.0)

    def test_censored_points(self, gal1_out: str) -> None:
        r = parse_km(gal1_out)
        assert len(r.censored_points) == 3
        assert 27.6 in r.censored_points
        assert 28.1 in r.censored_points
        assert 29.7 in r.censored_points

    def test_no_percentiles_when_too_few(self, gal1_out: str) -> None:
        r = parse_km(gal1_out)
        assert r.percentiles is None

    def test_differential(self, gal1_out: str) -> None:
        r = parse_km(gal1_out)
        assert r.differential is not None
        assert "bin_center" in r.differential.columns
        assert "D" in r.differential.columns
        assert len(r.differential) == 5

    def test_km_with_percentiles(self, gal2_out: str) -> None:
        # gal2.out contains two KM sections; the second (Starburst) has percentiles
        from asurv._parsers import parse_km
        import re

        km_sections = list(re.finditer(r"KAPLAN-MEIER ESTIMATOR", gal2_out))
        assert len(km_sections) >= 2
        second = gal2_out[km_sections[1].start():]
        r = parse_km(second)
        assert r.percentiles is not None
        assert 50 in r.percentiles
        assert abs(r.mean - 29.460) < 0.001


class TestParseTwoSample:
    def test_basic_fields(self, gal2_out: str) -> None:
        r = parse_twosample(gal2_out)
        assert r.n == 12
        assert r.n_upper_limits == 5
        assert r.n_group1 == 6
        assert r.n_group2 == 6

    def test_five_tests_present(self, gal2_out: str) -> None:
        r = parse_twosample(gal2_out)
        for key in ("gehan_perm", "gehan_hyper", "logrank", "peto_peto", "peto_prentice"):
            assert key in r.tests, f"missing test: {key}"

    def test_gehan_perm(self, gal2_out: str) -> None:
        r = parse_twosample(gal2_out)
        ts = r.tests["gehan_perm"]
        assert abs(ts.statistic - 1.652) < 0.001
        assert abs(ts.probability - 0.0986) < 0.0001

    def test_logrank(self, gal2_out: str) -> None:
        r = parse_twosample(gal2_out)
        ts = r.tests["logrank"]
        assert abs(ts.statistic - 1.814) < 0.001
        assert abs(ts.probability - 0.0696) < 0.0001

    def test_peto_prentice(self, gal2_out: str) -> None:
        r = parse_twosample(gal2_out)
        ts = r.tests["peto_prentice"]
        assert abs(ts.statistic - 1.706) < 0.001
        assert abs(ts.probability - 0.0881) < 0.0001

    def test_group_km_parsed(self, gal2_out: str) -> None:
        r = parse_twosample(gal2_out)
        assert r.km_group1 is not None
        assert r.km_group2 is not None
        assert r.km_group2.n == 6
        assert abs(r.km_group2.mean - 29.460) < 0.001


class TestParseBivar:
    def test_basic_fields(self, gal3_out: str) -> None:
        r = parse_bivar(gal3_out)
        assert r.n == 16
        assert r.upper_limits["y"] == 6
        assert r.upper_limits["x"] == 0

    def test_cox(self, gal3_out: str) -> None:
        r = parse_bivar(gal3_out)
        assert r.cox is not None
        assert abs(r.cox.chi2 - 18.458) < 0.001
        assert r.cox.dof == 1
        assert r.cox.probability < 0.001

    def test_em(self, gal3_out: str) -> None:
        r = parse_bivar(gal3_out)
        assert r.em is not None
        assert abs(r.em.intercept - 0.1703) < 0.001
        assert abs(r.em.intercept_err - 2.2356) < 0.001
        assert len(r.em.slopes) == 1
        slope, slope_err = r.em.slopes[0]
        assert abs(slope - 1.0519) < 0.001
        assert abs(slope_err - 0.0793) < 0.001
        assert abs(r.em.sigma - 0.3687) < 0.001
        assert r.em.iterations == 27

    def test_no_bhk_in_gal3(self, gal3_out: str) -> None:
        r = parse_bivar(gal3_out)
        assert r.bhk is None

    def test_no_bj_in_gal3(self, gal3_out: str) -> None:
        r = parse_bivar(gal3_out)
        assert r.bj is None
