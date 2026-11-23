#!/usr/bin/env python3
"""Unit test for ffz_casimir_lesage."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_casimir_lesage import (  # noqa: E402
    confinement_as_casimir,
    self_checks,
    verdict,
)


def test_self_checks():
    self_checks()


def test_avoids_lesage_not_two_body_theorem():
    v = verdict()
    assert v["avoids_all_lesage_fatals"] is True
    assert v["derived_two_body_Casimir"] is False
    assert confinement_as_casimir()["parallel_to_casimir_mode_subtraction"] is True


if __name__ == "__main__":
    test_self_checks()
    test_avoids_lesage_not_two_body_theorem()
    print("test_ffz_casimir_lesage: ALL PASS")
