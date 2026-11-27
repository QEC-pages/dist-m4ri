#!/usr/bin/env python3
"""Find every internal #anchor in every tracked .md file that pandoc fails to resolve."""
from __future__ import annotations

import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def gfm_slugs(lines: list[str]) -> dict[str, int]:
    res, seen, fence = {}, Counter(), False
    for idx, ln in enumerate(lines, 1):
        if ln.startswith("```"):
            fence = not fence
            continue
        m = None if fence else re.match(r"#{1,6} (.*)", ln)
        if m:
            s = re.sub(r"[^\w\- ]", "", m.group(1).strip().lower()).replace(" ", "-")
            slug = s if not seen[s] else f"{s}-{seen[s]}"
            seen[s] += 1
            res[slug] = idx
    return res


def main() -> int:
    files = subprocess.check_output(["git", "-C", str(ROOT), "ls-files", "*.md"], text=True).splitlines()
    for rel in files:
        p = ROOT / rel
        lines = p.read_text().splitlines()
        gmap = gfm_slugs(lines)
        tex = subprocess.check_output(["pandoc", str(p), "-t", "latex"], text=True)
        labels = set(re.findall(r"\\label\{([^}]+)\}", tex))
        hrefs = sorted(set(re.findall(r"\\hyperref\[([^\]]+)\]", tex)))
        missing = [h for h in hrefs if h not in labels]
        if missing:
            print(f"=== {rel}: {len(missing)} unresolved pandoc hyperrefs ===")
            for m in missing:
                lno = gmap.get(m)
                hdr = lines[lno - 1] if lno else "NOT EVEN IN GFM SLUGS!"
                print(f"  #{m} -> line {lno}: {hdr}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
