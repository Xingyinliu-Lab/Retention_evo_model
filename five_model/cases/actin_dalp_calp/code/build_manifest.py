#!/usr/bin/env python3
"""Create the final package manifest without modifying analytical outputs."""

import csv
import hashlib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
OUTPUT = PACKAGE / "MANIFEST_SHA256.tsv"
EXCLUDED_NAMES = {OUTPUT.name, "PACKAGE_COMPLETE.marker"}


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    rows = []
    for path in sorted(PACKAGE.rglob("*")):
        if not path.is_file() or path.name in EXCLUDED_NAMES:
            continue
        if "__pycache__" in path.parts:
            continue
        relative = path.relative_to(PACKAGE).as_posix()
        rows.append([relative, path.stat().st_size, digest(path)])
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["relative_path", "bytes", "sha256"])
        writer.writerows(rows)
    print(f"status=complete files={len(rows)}")


if __name__ == "__main__":
    main()
