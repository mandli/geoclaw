#!/usr/bin/env python
r"""Boussinesq (SGN) radial_flat test, run with PETSc under MPI.

Requires PETSc and MPI.  Skipped unless the environment variables used by the
example Makefile are set (see ``$CLAW/geoclaw/examples/bouss/setenv.sh``) and
point at an actual PETSc install and MPI commands.

Run with ``pytest -m bouss``.  This is not yet a regression comparison against
saved output; it checks that the run completes under MPI, that the MPI linear
solver server is used, and that water volume is conserved.
"""

import os
import re
import shutil
from pathlib import Path

import numpy as np
import pytest

import clawpack.geoclaw.test as test
import clawpack.pyclaw.solution as solution

_BOUSS_ENV = ("PETSC_DIR", "PETSC_OPTIONS", "CLAW_MPIEXEC", "CLAW_MPIFC")


def _bouss_unavailable() -> str:
    r"""Reason the Boussinesq build cannot run here, or "" if it can.

    Checks that PETSc and MPI are actually installed, not just that the
    environment variables are set (they may be left over in a shell profile).
    """
    missing = [name for name in _BOUSS_ENV if not os.environ.get(name)]
    if missing:
        return "PETSc/MPI environment not set: " + ", ".join(missing)
    petsc_dir = Path(os.environ["PETSC_DIR"])
    header = Path("petsc", "finclude", "petscksp.h")
    include_dirs = [petsc_dir / "include",
                    petsc_dir / os.environ.get("PETSC_ARCH", "") / "include"]
    if not any((d / header).is_file() for d in include_dirs):
        return f"PETSc Fortran headers not found under PETSC_DIR={petsc_dir}"
    for name in ("CLAW_MPIEXEC", "CLAW_MPIFC"):
        command = os.environ[name].split()[0]
        if shutil.which(command) is None:
            return f"{name}={os.environ[name]} not found on PATH"
    return ""


_UNAVAILABLE = _bouss_unavailable()
# Set in the Boussinesq CI job so a broken PETSc/MPI setup fails rather than
# silently skipping.
_REQUIRED = bool(os.environ.get("GEOCLAW_REQUIRE_BOUSS"))

pytestmark = [
    pytest.mark.regression,
    pytest.mark.bouss,
    pytest.mark.skipif(bool(_UNAVAILABLE) and not _REQUIRED,
                       reason=_UNAVAILABLE or "PETSc/MPI available"),
]


def _level1_stats(path: Path, frame: int) -> tuple[float, float]:
    r"""Total water volume and maximum depth on level 1."""
    sol = solution.Solution(frame, path=path, file_format="binary")
    volume = 0.0
    h_max = -np.inf
    for state in sol.states:
        if state.patch.level != 1:
            continue
        dx, dy = state.patch.delta
        volume += state.q[0, :, :].sum() * dx * dy
        h_max = max(h_max, state.q[0, :, :].max())
    return volume, h_max


def test_radial_flat(tmp_path: Path, monkeypatch):
    """Build and run the 2D Boussinesq radial_flat example under MPI."""

    if _UNAVAILABLE:
        pytest.fail(f"GEOCLAW_REQUIRE_BOUSS is set but {_UNAVAILABLE}")

    test_path = Path(__file__).parent
    num_procs = int(os.environ.get("BOUSS_MPI_PROCS", "2"))

    # setrun.py writes its flagregion files (RuledRectangle_*.data) into the
    # current directory and records their absolute paths, so run it from the
    # temporary directory.
    monkeypatch.chdir(tmp_path)

    runner = test.GeoClawTestRunner(tmp_path, test_path=test_path)
    runner.set_data()

    # Short run: the wave travels ~300 m of the 5000 m quarter domain, so
    # nothing reaches the extrapolation boundaries and the walls at x=0, y=0
    # carry no normal flux.
    runner.rundata.clawdata.tfinal = 10.0
    runner.rundata.clawdata.num_output_times = 1
    runner.rundata.clawdata.checkpt_style = 0
    runner.rundata.topo_data.topofiles = [[1, str(test_path / "flat100.tt1")]]
    runner.write_data()

    # The example Makefile sets EXE = $(PWD)/xgeoclaw, and PWD comes from the
    # environment rather than the directory make runs in; pin it.
    runner.build_executable(
        make_vars={"EXE": str(test_path / runner.executable_name)})
    runner.run_code(
        runexe=f"{os.environ['CLAW_MPIEXEC']} -n {num_procs}")

    # Mass conservation: the Boussinesq correction modifies only the momenta,
    # and with no boundary flux the total volume on level 1 (which holds the
    # conservative average of finer levels) should be unchanged to rounding.
    v0, h_max0 = _level1_stats(tmp_path, 0)
    v1, h_max1 = _level1_stats(tmp_path, 1)
    rel_change = abs(v1 - v0) / v0
    print(f"level-1 volume {v0!r} -> {v1!r} (rel {rel_change:.3e}); "
          f"max h {h_max0:.4f} -> {h_max1:.4f}")

    # Guard against a vacuous pass: the 5 m hump on 100 m of water must have
    # spread out, so the run actually advanced.
    assert h_max1 < h_max0 - 1.0, \
        f"solution barely changed: max h {h_max0} -> {h_max1}"
    assert rel_change <= 1e-12, \
        f"level-1 volume changed: {v0!r} -> {v1!r} (rel {rel_change:.3e})"

    # amr2.f90 redirects PETSc's stdout to timing.txt.  Check the MPI linear
    # solver server actually distributed at least one solve over all ranks,
    # so this test exercises the MPI path rather than a sequential fallback.
    timing = (tmp_path / "timing.txt").read_text()
    assert "MPI linear solver server statistics" in timing
    if num_procs > 1:
        assert re.search(rf"^\s*{num_procs}\s+[1-9]\d*\s", timing, re.MULTILINE), \
            f"no KSPSolve()s on {num_procs} ranks in timing.txt"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
