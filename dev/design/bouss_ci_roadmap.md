# Boussinesq (PETSc) CI — Roadmap

*Status: living document. Plan for testing the 2D Boussinesq solvers
(`src/2d/bouss`, PETSc + MPI) in GitHub Actions.*

**Status (2026-10-16): Phase 0 done; Phase 1 open as clawpack/clawutil#210;
Phase 2 implemented on `bouss-petsc-ci` and passing in a local replica of the
CI container (not yet run on Actions); Phase 3 not started; Phase 4 waiting on a
test case from the bouss authors (clawpack/geoclaw#759).**

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

**Implemented** on clawutil branch `test-runner-runexe` (base `027798c`):
`run_code(self, runexe: Optional[str]=None)` passes `runexe` to `runclaw`;
`None` keeps the old behavior. Because `run_example_for_test` already forwards
`run_kwargs` to `run_code`, a test can use either
`runner.run_code(runexe="mpiexec -n 2")` or
`run_example_for_test(..., run_kwargs={"runexe": ...})`.

Checked end to end through `runclaw`'s real `shlex`/`subprocess` path, using a
stand-in executable that records an environment variable set only by the
prefix (`/usr/bin/env PROBE=1`): the default run sees no prefix; `runexe` is
applied and runs in `temp_path`; the `run_kwargs` route forwards it. No
existing caller passes arguments to `run_code`, and no runner subclass
overrides it.

**Landing order (verified 2026-10-16).** CI checks out `clawpack/clawpack@master`
*with submodules*, which pins clawutil to the superrepo's recorded commit
(currently `027798cd`), not clawutil's master. So Phase 2 needs both:

1. the clawutil PR merged, **then**
2. a submodule bump of `clawutil` in `clawpack/clawpack` (precedent: "Bump
   clawutil for testing", `632541bc`).

## 3. Phase 2 — CI job and marker

**Implemented** on `bouss-petsc-ci`:

- **Marker:** `bouss` registered in `pytest.ini`.
- **Skip logic** (`examples/bouss/radial_flat/test_radial_flat.py`): skips
  unless the Makefile's variables are set **and** PETSc's Fortran header
  (`petsc/finclude/petscksp.h`) and the MPI commands actually exist. Env vars
  alone were not enough: a stale shell profile can set `PETSC_DIR` to a missing
  install, which made `make new` fail instead of skipping.
- **No exclusion needed in existing jobs** (changed from the original plan): the
  regression selectors collect the test, and it skips there because those jobs
  have no PETSc.
- **Fail, don't skip, in the bouss job:** the job sets `GEOCLAW_REQUIRE_BOUSS=1`,
  which turns the skip into a failure, so a broken PETSc setup can't go green.
- **Job** (`.github/workflows/bouss-petsc.yml`, replaces the probe): runs
  `-m bouss` in `ubuntu:26.04` with a `[debug, opt]` matrix, using
  `set_test_env.sh` for `FFLAGS` (the Makefile's `FFLAGS ?=` then
  `+= -DHAVE_PETSC` composes with them). `BOUSS_MPI_PROCS=2`. Uploads
  `pytest_tmp` on failure.
- **Triggers:** push/PR to `master` filtered to `src/2d/bouss/**`,
  `examples/bouss/**` and the workflow file; nightly at 06:30 UTC (catches the
  shared `shallow/`/amrclaw sources the filter can't see); `workflow_dispatch`.
- **Temporary:** a step checks out `mandli/clawutil@test-runner-runexe` over the
  pinned submodule. Remove it once clawutil#210 merges and the submodule is
  bumped.
- **Makefile trap:** `radial_flat/Makefile` sets `EXE = $(PWD)/xgeoclaw`, and
  `PWD` comes from the caller's environment, not make's working directory. The
  test pins `EXE` via the runner's `make_vars`.

**What the test checks** (not yet a comparison against saved output):

1. the build and an `mpiexec -n 2` run complete;
2. `timing.txt` shows the MPI linear solver server ran solves on all ranks;
3. **volume conservation** on level 1 between t=0 and t=10 (walls at x=0 and
   y=0, wave ~300 m from the far boundaries, and the Boussinesq correction only
   changes momenta), to a relative 1e-12;
4. a guard against a vacuous pass: max depth must fall by more than 1 m.

**Verified** in a local replica of the job (`ubuntu:26.04` arm64, clawutil from
`test-runner-runexe`), both builds with a fresh `make new` of all 162 sources:
volume 2500039269.9076505 → 2500039269.907651 (relative 1.9e-16), max depth
104.27 → 100.62 m, level-1 solves on 2 ranks. Debug and opt results were
bit-identical. About 6–7 s per build+run.

## 4. Phase 3 — saved-output regression comparison

Interim test from `examples/bouss/radial_flat`, cut down (few steps, one or two
output frames, a gauge or two), regression data generated once with `--save`.

- **Tolerance:** GMRES + GAMG at `ksp_rtol 1e-9` is not bit-reproducible across
  PETSc versions or rank counts. Pin ranks; compare with a relative tolerance
  (order 1e-8; to be set from observed run-to-run spread, not guessed).
- This is a **consistency** check only: it detects change, not correctness.

## 5. Phase 4 — a real test case (bouss authors)

Requested from Randy LeVeque and Marsha Berger in clawpack/geoclaw#759 (opened 2026-10-16).
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
- ~~Later, optional: macOS (Homebrew `petsc`) nightly cell.~~ Dropped
  2026-10-16: Homebrew `petsc` (3.26) requires `hdf5-mpi`, which conflicts
  with the `hdf5` that `netcdf`/`netcdf-fortran` build on, so it can't be
  installed alongside a normal GeoClaw NetCDF setup. On a Mac, use the CI
  container for bouss work.
- `radial_flat/Makefile`: `EXE = $(PWD)/xgeoclaw` → `$(CURDIR)/xgeoclaw`
  (`PWD` is the caller's, not make's, directory).
- **HPC (e.g. Derecho) — unverified, notes for users and docs:**
  - PETSc must be built with the *same compiler and MPI* as GeoClaw (Fortran
    `.mod` compatibility, as in CI). Use a site module only if it matches the
    loaded compiler/MPI and is ≥ 3.23; otherwise build PETSc in `$WORK`
    against the module toolchain. `petscMPIoptions` uses GMRES + GAMG, both
    built in, so no external packages are needed beyond BLAS/LAPACK.
  - The example Makefile's flags are gfortran-only (`-ffree-line-length-none`
    is appended unconditionally; the default `FFLAGS` has `-ffpe-trap`,
    `-gno-strict-dwarf`). Builds with Intel `ifx` are untested; the
    Makefile should get compiler-conditional flags.
  - Launching: rank 0 runs the AMR code (with OpenMP threads) and the other
    ranks only serve the linear solves (`-mpi_linear_solver_server`), so
    `BOUSS_MPI_PROCS` and `OMP_NUM_THREADS` must be placed together in the
    PBS allocation; single-node use is the intended case.

## 7. Risks

- **apt PETSc drift:** a 26.04 point update can change iteration counts and
  low-order digits. Tolerances must absorb it; the diagnostics step records the
  version in every log.
- **Image tag:** `ubuntu:26.04` moves with security updates; pin a digest only
  if drift becomes a problem.
- **Cross-repo ordering:** Phase 1 (clawutil) gates Phase 2.
