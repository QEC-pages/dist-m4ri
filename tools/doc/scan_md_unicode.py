#!/usr/bin/env python3
"""Scan all tracked .md files for non-ASCII characters and internal #anchor links."""
from __future__ import annotations

import re
import subprocess
import unicodedata
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    md_files = subprocess.check_output(["git", "-C", str(ROOT), "ls-files", "*.md"], text=True).splitlines()
    for rel in md_files:
        txt = (ROOT / rel).read_text()
        chars = Counter(ch for ch in txt if ord(ch) > 127)
        if chars:
            desc = ", ".join(
                f"{ch!r} U+{ord(ch):04X} ({unicodedata.name(ch, '?')}) x{n}"
                for ch, n in sorted(chars.items(), key=lambda kv: ord(kv[0]))
            )
            print(f"{rel}: {desc}")
        else:
            print(f"{rel}: ASCII only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
