#!/usr/bin/env python3
"""Add <div id="<gfm-slug>"></div> anchors above headings whose GFM slug is used by intra-file links
and differs from Pandoc's auto_identifiers label."""
from __future__ import annotations

import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def gfm_slugs(lines: list[str]) -> dict[str, int]:
    res, seen, fence = {}, Counter(), False
    for idx, ln in enumerate(lines):
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
        if not missing:
            continue
        inserts = sorted(((gmap[m], m) for m in missing), reverse=True)
        for idx, slug in inserts:
            tag = f'<div id="{slug}"></div>'
            # Ensure blank line before and after <div id="..."></div> so both GFM and Pandoc parse ## cleanly
            prefix = [tag, ""] if idx == 0 or lines[idx - 1] == "" else ["", tag, ""]
            lines[idx:idx] = prefix
        p.write_text("\n".join(lines) + "\n")
        print(f"Updated {rel}: inserted {len(inserts)} <div id=...> anchors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
