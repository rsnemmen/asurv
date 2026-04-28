"""Executable discovery, workspace lifecycle, and subprocess invocation."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


_KEEP_WORKSPACE_ENV = "ASURV_KEEP_WORKSPACE"
_BIN_ENV = "ASURV_BIN"


def find_executable() -> Path:
    """Return the path to the asurv binary.

    Search order:
    1. ASURV_BIN environment variable
    2. ./asurv relative to the current working directory
    3. asurv on PATH
    """
    from asurv import AsurvNotFoundError

    if env := os.environ.get(_BIN_ENV):
        p = Path(env)
        if p.is_file():
            return p
        raise AsurvNotFoundError(
            f"ASURV_BIN={env!r} is set but the file does not exist"
        )

    local = Path.cwd() / "asurv"
    if local.is_file():
        return local

    on_path = shutil.which("asurv")
    if on_path:
        return Path(on_path)

    raise AsurvNotFoundError(
        "asurv executable not found. "
        "Run `make` in the repo root, then retry from that directory, "
        "or set the ASURV_BIN environment variable to its path."
    )


class Workspace:
    """Temporary directory that is cleaned up unless ASURV_KEEP_WORKSPACE=1."""

    def __init__(self) -> None:
        self.path = Path(tempfile.mkdtemp(prefix="asurv_"))

    def __enter__(self) -> "Workspace":
        return self

    def __exit__(self, *_: object) -> None:
        keep = os.environ.get(_KEEP_WORKSPACE_ENV, "").strip() not in ("", "0")
        if keep:
            import sys
            print(f"asurv workspace kept at {self.path}", file=sys.stderr)
        else:
            shutil.rmtree(self.path, ignore_errors=True)


def run(
    mode: str,
    data_text: str,
    command_lines: list[str],
    data_filename: str = "data.dat",
    command_filename: str | None = None,
    output_filename: str | None = None,
) -> str:
    """Invoke the asurv executable in batch mode and return the output text.

    Parameters
    ----------
    mode:
        One of ``"KM"``, ``"TWOST"``, or ``"BIVAR"``.
    data_text:
        Contents of the data file.
    command_lines:
        Lines of the command file (without a trailing newline on each).
    data_filename:
        Name written inside the workspace (≤9 chars, kept short by callers).
    command_filename:
        Command file name inside the workspace (defaults to ``<mode.lower()>.com``).
    output_filename:
        Expected output file name (defaults to ``<mode.lower()>.out``).
        Pass ``None`` to use the default.
    """
    from asurv import AsurvExecutionError

    if command_filename is None:
        command_filename = f"{mode.lower()}.com"
    if output_filename is None:
        output_filename = f"{mode.lower()}.out"

    exe = find_executable()

    with Workspace() as ws:
        (ws.path / data_filename).write_text(data_text)
        (ws.path / command_filename).write_text("\n".join(command_lines) + "\n")

        completed = subprocess.run(
            [str(exe.resolve()), mode, command_filename],
            cwd=ws.path,
            capture_output=True,
            text=True,
            check=False,
        )

        output_path = ws.path / output_filename
        if completed.returncode != 0:
            raise AsurvExecutionError(
                f"asurv exited with status {completed.returncode}:\n"
                f"{completed.stdout}{completed.stderr}"
            )
        if not output_path.is_file():
            raise AsurvExecutionError(
                f"asurv did not produce expected output file {output_filename!r}.\n"
                f"{completed.stdout}{completed.stderr}"
            )
        return output_path.read_text()
