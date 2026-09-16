"""Central numerical thresholds from the frozen policy."""

import json
from functools import lru_cache

from .config import repo_root


@lru_cache(maxsize=1)
def numerical_policy():
    return json.loads((repo_root() / "NUMERICAL_POLICY.json").read_text())
