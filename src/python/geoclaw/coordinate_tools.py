#!/usr/bin/env python
# encoding: utf-8
r"""
Format- and product-neutral coordinate helpers for GeoClaw input.

This module holds the longitude/geometry primitives that are shared by *all*
input products (topography, moving topography / dtopo, and gridded met forcing)
and *all* file formats (ASCII types 1/2/3 and NetCDF type 4).  Keeping them here
-- rather than in ``netcdf_utils`` -- gives every path a single implementation
to import, so the antimeridian-wrap logic is never copied per format.

Contents
--------
is_geographic_lon
    Decide whether an x axis is a geographic longitude (so the +/-360 wrap is
    meaningful) or a projected / Cartesian axis (which never wraps).
_compute_lon_entries
    Pure geometry: given a file's longitude span and a requested domain span,
    return the ``(file_crop_min, file_crop_max, lon_offset)`` tuples needed to
    cover the domain -- one tuple normally, two when the request straddles the
    file's longitude cut (dateline).
"""

from __future__ import annotations

from typing import Optional

import numpy


# Length unit strings (lower-cased) on an x axis that mark it as a projected /
# rectilinear (non-geographic) grid, for which the 0-360 longitude wrap is
# meaningless and must be skipped.  Detection defaults to geographic; only
# positive evidence (these units, a projection standard_name, or values
# outside the plausible degree range) disables the wrap.
_PROJECTED_LENGTH_UNITS: frozenset[str] = frozenset({
    'm', 'meter', 'meters', 'metre', 'metres',
    'km', 'kilometer', 'kilometers', 'kilometre', 'kilometres',
})
_PROJECTED_STANDARD_NAMES: frozenset[str] = frozenset({
    'projection_x_coordinate', 'projection_y_coordinate',
})


def is_geographic_lon(
    lon_min: float,
    lon_max: float,
    units: Optional[str] = None,
    std_name: Optional[str] = None,
) -> bool:
    """Return True if an x axis should be treated as geographic longitude.

    The axis is assumed geographic unless there is positive evidence it is a
    projected / Cartesian (e.g. beta-plane) grid: a length unit (``m``/``km``),
    a projection ``standard_name``, or coordinate values outside the plausible
    degree range ``[-360, 360]``.  Antimeridian wrapping (+/-360) is only ever
    applied when this returns True; a Cartesian grid is left untouched.

    ``units``/``std_name`` are optional (ASCII files carry no such metadata);
    when omitted, only the coordinate-magnitude test applies.
    """
    u = str(units or '').strip().lower()
    s = str(std_name or '')
    projected = (
        u in _PROJECTED_LENGTH_UNITS
        or s in _PROJECTED_STANDARD_NAMES
        or lon_min < -360.0 - 1e-6
        or lon_max > 360.0 + 1e-6
    )
    return not projected


def _compute_lon_entries(
    file_lon_min: float,
    file_lon_max: float,
    domain_lon_min: float,
    domain_lon_max: float,
    max_gap: float = 1e-10,
    allow_wrap: bool = True,
) -> list[tuple[float, float, float]]:
    """
    Compute (file_crop_min, file_crop_max, lon_offset) tuples needed to cover
    [domain_lon_min, domain_lon_max] from a file with lons in
    [file_lon_min, file_lon_max].

    lon_offset is the scalar added to file coordinates to produce domain
    coordinates: x_domain = x_file + lon_offset.  For NetCDF this is the
    descriptor ``lon_wrap_offset``; for ASCII it is folded into ``x_shift``.

    Returns 1 tuple if a single offset suffices, 2 tuples if the domain
    straddles the file's cut point.

    max_gap controls how much under-coverage is tolerated.  The default
    (1e-10) is tight enough to catch genuine gaps.  Pass the file's grid
    spacing to allow for the half-cell gap at the dateline that near-global
    files (e.g. GEBCO) have between their last and first longitude columns.

    allow_wrap enables the +/-360 wrap candidate offsets (the default, for a
    geographic longitude axis).  Pass False for a non-geographic x axis
    (projected meters, etc.), which never wraps: only the identity offset 0
    is considered.

    Raises ValueError if the file cannot cover the requested domain even
    with all candidate offsets, or if the remaining uncovered gap exceeds
    max_gap.
    """
    candidate_offsets = [0.0, 360.0, -360.0] if allow_wrap else [0.0]
    entries: list[tuple[float, float, float]] = []
    total_coverage = 0.0

    for offset in candidate_offsets:
        shifted_min = file_lon_min + offset
        shifted_max = file_lon_max + offset
        intersect_min = max(domain_lon_min, shifted_min)
        intersect_max = min(domain_lon_max, shifted_max)
        width = intersect_max - intersect_min
        if width > 1e-10:
            file_crop_min = intersect_min - offset
            file_crop_max = intersect_max - offset
            # Clamp to file extent (guards against floating-point overshoot)
            file_crop_min = max(file_crop_min, file_lon_min)
            file_crop_max = min(file_crop_max, file_lon_max)
            entries.append((file_crop_min, file_crop_max, offset))
            total_coverage += width

    if not entries:
        raise ValueError(
            f"File longitude range [{file_lon_min}, {file_lon_max}] cannot cover "
            f"domain [{domain_lon_min}, {domain_lon_max}] with candidate offsets "
            f"{candidate_offsets}."
        )

    # One-sided check: overcoverage is harmless; under-coverage beyond max_gap
    # indicates that the file genuinely cannot cover the requested domain.
    gap = (domain_lon_max - domain_lon_min) - total_coverage
    if gap > max_gap:
        raise ValueError(
            f"File longitudes [{file_lon_min}, {file_lon_max}] cannot "
            f"cover requested domain [{domain_lon_min}, {domain_lon_max}]. "
            f"Gap of {gap:.6f} degrees exceeds tolerance {max_gap:.6f}."
        )

    return entries


def crop_indices(coords, lo, hi, delta, coarsen=1, buffer=0, align=None):
    r"""Half-open slice bounds ``(lower, upper)`` for a 1-D crop + coarsen.

    ``coords[lower:upper:coarsen]`` is the cropped, coarsened, buffered,
    align-adjusted lattice covering the inclusive request ``[lo, hi]``.  This is
    the single index-window implementation shared by the topography and moving-
    topography (dtopo) crop paths, for both ASCII and NetCDF products.

    :Input:
     - *coords* (ndarray) - 1-D ascending coordinate array with uniform spacing.
     - *lo*, *hi* (float) - inclusive requested crop bounds (domain coords).
     - *delta* (float) - grid spacing of *coords* (``coords[1] - coords[0]``).
     - *coarsen* (int) - coarsening (subsampling) factor; 1 = none.
     - *buffer* (int) - integer number of grid points to keep on each side (NOT
       a coordinate distance).  Scaled by *coarsen* so it counts points on the
       original, pre-coarsen lattice.
     - *align* (float or None) - alignment target for coarsening; when given and
       ``coarsen > 1`` the low index is shifted so ``(coords[lower] - align) /
       (delta*coarsen)`` is as close to an integer as the lattice allows.

    Returns ``(lower, upper)`` (Python ints, half-open) or ``None`` when
    ``[lo, hi]`` does not overlap *coords* -- callers preserve their own
    "does not overlap" handling.
    """
    coords = numpy.asarray(coords)
    n = len(coords)
    try:
        lower = (coords >= lo).nonzero()[0][0]
        upper = (coords <= hi).nonzero()[0][-1]
    except IndexError:
        return None

    delta_new = delta * coarsen

    # Shift the low index for alignment when coarsening (mirrors the original
    # Topography.crop math exactly: search the first `coarsen` points for the
    # best-aligned start, then trim the high index to a whole coarsen stride).
    if (coarsen > 1) and (align is not None):
        vs = numpy.array([coords[lower + i] for i in range(coarsen)])
        offsets = (vs - align) / delta_new
        offsets_frac = offsets - numpy.round(offsets)
        ioffset = numpy.argmin(abs(offsets_frac))
        lower = lower + ioffset
        upper = upper - numpy.remainder(upper - lower, coarsen)

    # Buffer (integer grid-point count), clamped to the array limits.
    lower = numpy.maximum(0, lower - buffer * coarsen)
    upper = numpy.minimum(n - 1, upper + buffer * coarsen) + 1

    return int(lower), int(upper)


def netcdf_window_indices(coords, lo, hi, margin, n, stride=1):
    r"""Half-open index window ``[i0, i1)`` into 1-D monotonic *coords*.

    A read-time *prefilter* (distinct from :func:`crop_indices`): it pushes a
    requested crop down to a NetCDF hyperslab read so only the needed slab is
    loaded from disk, keeping a *margin* of extra points on each side so a
    subsequent exact :func:`crop_indices` still has every point it needs.  It is
    format-neutral geometry, shared by every NetCDF input reader.

    :Input:
     - *coords* (ndarray) - 1-D coordinate array, monotonic ascending **or**
       descending (as stored in the file).
     - *lo*, *hi* (float) - requested inclusive coordinate bounds.
     - *margin* (int) - extra points to keep on each side (the caller's buffer
       plus the coarsen alignment search room).
     - *n* (int) - length of *coords*.
     - *stride* (int) - read stride; the low index is snapped **down** to a
       multiple of *stride* so the strided sub-window shares the same phase as
       striding the full array, keeping the sampled grid identical.

    Returns ``(0, n)`` (the full range) when the interval does not overlap
    *coords*, mirroring the crop fall-back of leaving the array uncropped when
    the requested region misses the data.
    """
    # A monotonic array clipped to [lo, hi] yields a contiguous True block for
    # either sort order, so first/last True index bound the window.
    mask = (coords >= lo) & (coords <= hi)
    idx = numpy.nonzero(mask)[0]
    if idx.size == 0:
        return 0, n
    i0 = max(0, int(idx[0]) - margin)
    i1 = min(n, int(idx[-1]) + margin + 1)
    # Snap the low index down to a multiple of stride so coords[i0:i1:stride]
    # is a phase-aligned subset of coords[::stride].
    i0 -= i0 % stride
    if i1 <= i0:
        return 0, n
    return i0, i1
