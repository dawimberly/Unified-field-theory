#!/usr/bin/env python3
"""Unit test for ffz_string_kk_tower."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_string_kk_tower import proca_tower_flip, self_checks, species_bound_note  # noqa: E402


def test_self_checks():
    self_checks()


def test_proca_count():
    assert proca_tower_flip()["N_Proca_min_to_flip_UV"] == 8
    assert species_bound_note()["guarantees_positive_a1"] is False


if __name__ == "__main__":
    test_self_checks()
    test_proca_count()
    print("test_ffz_string_kk_tower: ALL PASS")
