"""Portable paths for the migrated numerical scripts; archived originals are untouched."""

from .config import repo_root


class RepositoryPaths:
    def __truediv__(self, name):
        root = repo_root()
        for folder in (
            root / "src",
            root / "data/authoritative",
            root / "data/regression",
            root / "data/certificates",
            root / "data/production",
            root,
        ):
            path = folder / name
            if path.is_file():
                return path
        output = root / "build" / name
        output.parent.mkdir(parents=True, exist_ok=True)
        return output


paths = RepositoryPaths()
