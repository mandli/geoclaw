#!/usr/bin/env python
# encoding: utf-8
"""Tests for the product-/format-neutral coordinate primitives."""

import pytest

from clawpack.geoclaw.coordinate_tools import (
    is_geographic_lon,
    _compute_lon_entries,
)


@pytest.mark.python
def test_is_geographic_lon_degree_ranges():
    assert is_geographic_lon(-190.0, -130.0) is True   # extended-frame lon-lat
    assert is_geographic_lon(0.0, 360.0) is True
    assert is_geographic_lon(-180.0, 180.0) is True


@pytest.mark.python
def test_is_geographic_lon_projected():
    # Large (meters-like) coordinates -> projected, never wrap.
    assert is_geographic_lon(0.0, 500000.0) is False
    # Length units on the axis -> projected even for small coords.
    assert is_geographic_lon(0.0, 10.0, units="km") is False
    assert is_geographic_lon(0.0, 10.0, units="m") is False
    # Projection standard_name -> projected.
    assert is_geographic_lon(0.0, 10.0,
                             std_name="projection_x_coordinate") is False
    # Small unit-less coords default to geographic (heuristic can't tell).
    assert is_geographic_lon(0.0, 10.0) is True


@pytest.mark.python
def test_compute_lon_entries_single_and_split():
    # Non-crossing crop within a [-180,180] file -> one identity entry.
    entries = _compute_lon_entries(-180.0, 180.0, -170.0, -150.0)
    assert entries == [(-170.0, -150.0, 0.0)]

    # Crop crossing -180 -> two entries with offsets {0, -360}.
    entries = _compute_lon_entries(-180.0, 180.0, -190.0, -170.0, max_gap=1.0)
    offsets = sorted(o for *_, o in entries)
    assert offsets == [-360.0, 0.0]
    pairs = {(round(a, 6), round(b, 6), o) for a, b, o in entries}
    assert (-180.0, -170.0, 0.0) in pairs
    assert (170.0, 180.0, -360.0) in pairs


@pytest.mark.python
def test_compute_lon_entries_no_wrap_disables_split():
    # allow_wrap=False (projected axis): only the identity offset is tried, so a
    # crop reaching past the file cannot be covered -> ValueError.
    with pytest.raises(ValueError):
        _compute_lon_entries(-180.0, 180.0, -190.0, -170.0,
                             max_gap=1.0, allow_wrap=False)
