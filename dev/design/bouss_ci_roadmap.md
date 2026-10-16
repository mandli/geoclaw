# Boussinesq (PETSc) CI — Roadmap

*Status: living document. Plan for testing the 2D Boussinesq solvers
(`src/2d/bouss`, PETSc + MPI) in GitHub Actions.*

**Status (2026-10-16): Phase 0 done; Phases 1–2 not started; Phase 4 waiting on
a test case from the bouss authors.** The toolchain probe runs green on GitHub
Actions (see §1). No bouss *test* exists yet.

**Why.** The bouss code compiles a separate source list (`Makefile.bouss`) that
mixes its own copies of AMR routines with shared `amrclaw/src/2d` and
`geoclaw/src/2d/shallow` sources, so changes elsewhere can break it, and CI has
never built it. `tests/regression/force_dry/test_force_dry_index.py` notes that
the `src/2d/bouss/` copies of `filpatch.f90`/`filval.f90` cannot be exercised
because CI lacks PETSc and MPI.

---

## 1. Phase 0 — toolchain probe (done)

`.github/workflows/bouss-petsc.yml` on branch `bouss-petsc-ci` (push-to-branch
and `workflow_dispatch` only). Builds `examples/bouss/radial_flat` and does a
short (`tfinal=10`) run under `mpiexec -n 2`, with no comparison against
reference output.

**Verified** against commit `59148c03`, mandli/geoclaw Actions run
38910482227 (all steps green, job 2m47s); also run locally in Docker on
arm64 and emulated amd64 with identical results.

Facts established (each observed, not assumed):

| Item | Finding |
|---|---|
| PETSc source | `ubuntu:26.04` container, apt `petsc-dev` → **3.24.4**. Ubuntu 24.04 ships 3.19, 25.10 ships 3.22; the requirement is ≥ 3.23 (`examples/bouss/README.txt`). |
| Compiler | Must be the image's gfortran (15.2). PETSc's `.mod` files are gfortran module version 16; `setup-fortran`'s gcc 14 is not used. |
| `PETSC_DIR` / `PETSC_ARCH` | `/usr/lib/petsc`, `PETSC_ARCH` empty. The example Makefile's `pkg-config` path resolves to `-lpetsc_real` unchanged. |
| MPI | Open MPI. Container needs `OMPI_ALLOW_RUN_AS_ROOT=1`, `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1`, `--shm-size=1g`. |
| Solver server | Engaged: level 1 (20 000 rows) solved on 2 ranks, ~8 GMRES its; level 2 (~4 000 rows, below `minimum_count_per_rank 5000`) sequential. |
| Where PETSc output goes | `amr2.f90` calls `PetscViewerASCIIStdoutSetFileUnit(timing_unit)`, so `-options_view` and solver-server stats land in **`timing.txt`**, not stdout. |
| Python deps | `setrun` write needs `xarray` even for a `.tt1` topo (`data.py` → `netcdf_utils`). |
| Container plumbing | Install `git` before `actions/checkout` (else no submodules); `${{ github.workspace }}` is the host path, so set `CLAW` from `$GITHUB_WORKSPACE`. |
| Makefile noise | `FC=mpif90` is not recognized by `Makefile.common` (warning, falls back to `-J`/`-fopenmp` — correct). |

## 2. Phase 1 — runner support (clawutil, lands first)

`ClawpackTestRunner.run_code` (`clawutil/src/python/clawutil/test.py`) calls
`runclaw.runclaw(...)` without `runexe`, so a test cannot launch under
`mpiexec`. `runclaw` already accepts `runexe`.

- Add an optional `runexe` (attribute or `run_code` argument, default `None`
  → current behavior) and pass it through.
- clawutil PR; must merge before Phase 2 because CI checks out
  `clawpack/clawpack@master`.

## 3. Phase 2 — CI job and marker

- **Marker:** add `bouss` (and/or `petsc`) to `pytest.ini`. Tests also
  `skipif` when `PETSC_DIR` is unset, so a plain local `pytest` doesn't fail.
- **Exclusion from existing jobs:** the existing regression jobs select
  `-m "regression and ..."`; add `and not bouss` there so they don't attempt it.
- **Job:** replace the probe with a job running `-m bouss` in the
  `ubuntu:26.04` container. Reuse `.github/scripts/set_test_env.sh` for
  `FFLAGS` (the example Makefile uses `FFLAGS ?=` then `+= -DHAVE_PETSC`, so
  exported flags compose). Pin `BOUSS_MPI_PROCS=2`, `OMP_NUM_THREADS`.
- **Triggers:** nightly (in `slow-tests.yml` or its own schedule) **plus**
  per-PR with a `paths:` filter on `src/2d/bouss/**`, `examples/bouss/**` and
  the workflow file. The filter cannot cover shared `shallow/` / amrclaw
  sources that `Makefile.bouss` pulls in; nightly is what catches those.
- On failure upload `fort.*`, `timing.*`, `*.data`.

## 4. Phase 3 — smoke regression test

Interim test from `examples/bouss/radial_flat`, cut down (few steps, one or two
output frames, a gauge or two), regression data generated once with `--save`.

- **Tolerance:** GMRES + GAMG at `ksp_rtol 1e-9` is not bit-reproducible across
  PETSc versions or rank counts. Pin ranks; compare with a relative tolerance
  (order 1e-8; to be set from observed run-to-run spread, not guessed).
- This is a **consistency** check only: it detects change, not correctness.

## 5. Phase 4 — a real test case (bouss authors)

Requested from Randy LeVeque and Marsha Berger (draft request in `.plans/`).
What we want it to exercise, in priority order:

1. Something with a known answer: 2D `radial_flat` vs. the `1d_radial`
   reference solution already in the example, or a linear-dispersion check
   against the SGN (`alpha = 1.153`) dispersion relation.
2. AMR coupling of the implicit solve: a fine patch boundary inside the
   Boussinesq region.
3. The SGN→SWE switch (`bouss_min_depth`) with topography crossing it.
4. Restart (bouss copies of `check.f` / `restrt.f`).
5. Optionally `bouss_min_level > 1`.

Constraint: total run (excluding build) well under a minute on a 4-core
runner.

## 6. Phase 5 — docs and follow-ups

- `clawpack/doc` `bouss2d.rst` still says PETSc 3.20; README says 3.23. Align.
- Document the Ubuntu 26.04 apt route (`petsc-dev`, `PETSC_DIR=/usr/lib/petsc`,
  no `PETSC_ARCH`) as a user install option.
- Update `test_force_dry_index.py`'s note once the bouss copies are covered.
- Later, optional: macOS (Homebrew `petsc`) nightly cell.

## 7. Risks

- **apt PETSc drift:** a 26.04 point update can change iteration counts and
  low-order digits. Tolerances must absorb it; the diagnostics step records the
  version in every log.
- **Image tag:** `ubuntu:26.04` moves with security updates; pin a digest only
  if drift becomes a problem.
- **Cross-repo ordering:** Phase 1 (clawutil) gates Phase 2.
