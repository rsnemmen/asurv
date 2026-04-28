"""Extract gal1/gal2/gal3 fixtures from asurv.etc for the pytest suite."""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pytest

_ASURV_ETC = Path(__file__).resolve().parents[1] / "asurv.etc"

_SECTION_RE = re.compile(
    r"\*{10,}\s*\*+\s*([\w.]+)\s*\*+\s*\*{10,}\n(.*?)(?=\*{10,}|\Z)",
    re.DOTALL,
)


def _parse_etc(path: Path) -> dict[str, str]:
    """Return a mapping of filename → content from asurv.etc."""
    text = path.read_text()
    sections: dict[str, str] = {}
    for m in _SECTION_RE.finditer(text):
        name = m.group(1).strip()
        content = m.group(2)
        # strip the leading/trailing blank lines added by the separator
        sections[name] = content
    return sections


@pytest.fixture(scope="session")
def asurv_fixtures(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Extract all gal*.dat/com/out files from asurv.etc into a session tmp dir."""
    sections = _parse_etc(_ASURV_ETC)
    base = tmp_path_factory.mktemp("fixtures")
    result: dict[str, Path] = {}
    for name, content in sections.items():
        if name.startswith("gal"):
            p = base / name
            p.write_text(content)
            result[name] = p
    return result


@pytest.fixture(scope="session")
def gal1_dat(asurv_fixtures: dict[str, Path]) -> Path:
    return asurv_fixtures["gal1.dat"]


@pytest.fixture(scope="session")
def gal2_dat(asurv_fixtures: dict[str, Path]) -> Path:
    return asurv_fixtures["gal2.dat"]


@pytest.fixture(scope="session")
def gal3_dat(asurv_fixtures: dict[str, Path]) -> Path:
    return asurv_fixtures["gal3.dat"]


@pytest.fixture(scope="session")
def gal1_out(asurv_fixtures: dict[str, Path]) -> str:
    return asurv_fixtures["gal1.out"].read_text()


@pytest.fixture(scope="session")
def gal2_out(asurv_fixtures: dict[str, Path]) -> str:
    return asurv_fixtures["gal2.out"].read_text()


@pytest.fixture(scope="session")
def gal3_out(asurv_fixtures: dict[str, Path]) -> str:
    return asurv_fixtures["gal3.out"].read_text()
