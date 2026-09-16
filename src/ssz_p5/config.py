from __future__ import annotations

from pathlib import Path

RELEASE = "2026-09-16"
EPSILON_Y = 1.0e-2
DEFAULT_L = (6, 12, 20, 42, 110, 420, 1000)
JET_WINDOW = 9
JET_DEGREE = 8
SLOT_NAMES = tuple(
    [f"a{i}" for i in range(1, 10)]
    + [f"b{i}" for i in range(1, 6)]
    + [f"c{i}" for i in range(1, 7)]
    + [f"d{i}" for i in range(1, 5)]
    + [f"e{i}" for i in range(1, 5)]
    + [f"v{i}" for i in range(1, 14)]
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]
