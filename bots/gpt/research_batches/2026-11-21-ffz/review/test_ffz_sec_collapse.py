#!/usr/bin/env python3
"""Unit test for ffz_sec_collapse."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_sec_collapse import (  # noqa: E402
    mit_bag_quark_matter,
    neutron_star_core,
    self_checks,
    verdict,
)


def test_self_checks():
    self_checks()


def test_only_vacuum_violates():
    assert neutron_star_core()["any_SEC_violation"] is False
    bag = mit_bag_quark_matter()
    assert bag["violation_only_bag_vacuum"] is True
    assert verdict()["ordinary_matter_violates_SEC"] is False


if __name__ == "__main__":
    test_self_checks()
    test_only_vacuum_violates()
    print("test_ffz_sec_collapse: ALL PASS")
