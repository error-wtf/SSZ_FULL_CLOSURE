"""Load migrated golden modules once; no numerical formulas are duplicated."""

import importlib.util
from functools import cache

from .config import repo_root


@cache
def module(filename):
    path = repo_root() / "src" / filename
    spec = importlib.util.spec_from_file_location(
        "ssz_p5_reference_" + filename.replace(".", "_"), path
    )
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj
