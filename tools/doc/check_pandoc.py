#!/usr/bin/env python3
"""Test compiling every tracked .md file with pandoc -> PDF and capture all errors/warnings."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    md_files = subprocess.check_output(["git", "-C", str(ROOT), "ls-files", "*.md"], text=True).splitlines()
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        any_issue = False
        for rel in md_files:
            p = ROOT / rel
            out_pdf = tmp_dir / (rel.replace("/", "__") + ".pdf")
            proc = subprocess.run(
                ["pandoc", p.name, "-o", str(out_pdf)],
                cwd=str(p.parent),
                capture_output=True,
                text=True,
            )
            stderr = proc.stderr.strip()
            status = "OK" if proc.returncode == 0 and not stderr else f"rc={proc.returncode}"
            print(f"=== {rel} [{status}] ===")
            if stderr:
                any_issue = True
                print(stderr)
            elif proc.returncode != 0:
                any_issue = True
        return 1 if any_issue else 0


if __name__ == "__main__":
    raise SystemExit(main())
