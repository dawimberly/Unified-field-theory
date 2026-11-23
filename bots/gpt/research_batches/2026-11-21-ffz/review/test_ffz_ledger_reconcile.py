#!/usr/bin/env python3
"""Unit test for ffz_ledger_reconcile."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_ledger_reconcile import channel_map, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_no_ledger_change():
    assert channel_map()["any_changes_ledger"] is False


if __name__ == "__main__":
    test_self_checks()
    test_no_ledger_change()
    print("test_ffz_ledger_reconcile: ALL PASS")
