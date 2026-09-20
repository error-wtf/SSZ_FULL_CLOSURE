#!/usr/bin/env python3
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def summarize(path: Path) -> tuple[list[str], list[str]]:
    try:
        tree = ast.parse(path.read_text())
    except (SyntaxError, UnicodeDecodeError):
        return [], []
    funcs = [n.name for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    classes = [n.name for n in tree.body if isinstance(n, ast.ClassDef)]
    return funcs, classes


def main() -> int:
    roots = [ROOT / "src" / "ssz_p5", ROOT / "tools"]
    for base in roots:
        print(f"\n## {base.relative_to(ROOT)}")
        for path in sorted(base.rglob("*.py")):
            funcs, classes = summarize(path)
            rel = path.relative_to(ROOT)
            bits = []
            if classes:
                bits.append("classes=" + ",".join(classes[:8]))
            if funcs:
                bits.append("functions=" + ",".join(funcs[:12]))
            print(f"{rel}: {'; '.join(bits) if bits else '(module)'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
