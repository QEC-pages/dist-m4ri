#!/usr/bin/env python3
"""Measure exact Overfull \\hbox warnings in pdflatex logs for all .md files."""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def find_longtable_ranges(tex_lines: list[str]) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    start = None
    for idx, line in enumerate(tex_lines, start=1):
        if r"\begin{longtable}" in line:
            start = idx
        elif r"\end{longtable}" in line and start is not None:
            ranges.append((start, idx))
            start = None
    return ranges


def in_ranges(line_no: int, ranges: list[tuple[int, int]]) -> bool:
    return any(lo <= line_no <= hi for lo, hi in ranges)


def main() -> int:
    geom = sys.argv[1] if len(sys.argv) > 1 else None
    verbose = "--verbose" in sys.argv
    if geom == "--verbose":
        geom = None
    md_files = subprocess.check_output(["git", "-C", str(ROOT), "ls-files", "*.md"], text=True).splitlines()
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for rel in md_files:
            p = ROOT / rel
            tex_file = tmp_dir / "doc.tex"
            cmd = ["pandoc", p.name, "-s", "-o", str(tex_file)]
            if geom:
                cmd += ["-V", f"geometry:{geom}"]
            subprocess.run(cmd, cwd=str(p.parent), check=True)
            subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-output-directory", str(tmp_dir), str(tex_file)],
                cwd=str(p.parent),
                capture_output=True,
                text=True,
                errors="replace",
            )
            log = (tmp_dir / "doc.log").read_text(errors="replace")
            tex_lines = tex_file.read_text().splitlines()
            lt_ranges = find_longtable_ranges(tex_lines)
            over = re.findall(
                r"Overfull \\hbox \(([0-9.]+)pt too wide\) in (alignment|paragraph) at lines ([0-9-]+)", log
            )
            tbl_over = []
            prose_over = []
            for pts_s, kind, lrange in over:
                pts = float(pts_s)
                lo = int(lrange.split("--")[0])
                if kind == "alignment" or in_ranges(lo, lt_ranges):
                    tbl_over.append((pts, kind, lo, lrange))
                elif pts > 5.0:
                    prose_over.append((pts, kind, lo, lrange))
            max_tbl = max((x[0] for x in tbl_over), default=0.0)
            max_prose = max((x[0] for x in prose_over), default=0.0)
            print(
                f"{rel:30s} table_over={len(tbl_over):3d} (max={max_tbl:6.1f}pt)  "
                f"prose_over(>5pt)={len(prose_over):3d} (max={max_prose:6.1f}pt)"
            )
            if verbose:
                for pts, kind, lo, lrange in tbl_over:
                    ctx = " ".join(tex_lines[max(0, lo - 2) : min(len(tex_lines), lo + 2)])[:110]
                    print(f"    [TBL {kind:9s}] +{pts:5.1f}pt @L{lrange}: {ctx}")
                for pts, kind, lo, lrange in prose_over[:10]:
                    ctx = " ".join(tex_lines[max(0, lo - 1) : min(len(tex_lines), lo + 2)])[:110]
                    print(f"    [PROSE]         +{pts:5.1f}pt @L{lrange}: {ctx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
