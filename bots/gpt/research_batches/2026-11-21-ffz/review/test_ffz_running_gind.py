#!/usr/bin/env python3
"""Unit test for ffz_running_gind."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_running_gind import crossover, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_crossover_at_qcd():
    assert crossover()["at_Lambda_QCD"] is True


if __name__ == "__main__":
    test_self_checks()
    test_crossover_at_qcd()
    print("test_ffz_running_gind: ALL PASS")
