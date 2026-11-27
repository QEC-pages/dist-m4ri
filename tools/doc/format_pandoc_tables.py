#!/usr/bin/env python3
"""Format wide pipe-table separator rows with proportional dash counts summing to 80."""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEP_RE = re.compile(r"^\|(\s*:?-+:?\s*\|)+$")


def split_row(line: str) -> list[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith(r"\|"):
        s = s[:-1]
    return [c.strip() for c in re.split(r"(?<!\\)\|", s)]


def approx_math_text(m: str) -> str:
    s = m
    s = re.sub(r"\\\|\\mathcal\{([A-Za-z0-9]+)\}_\{?([A-Za-z0-9]+)\}?\\\|", r"||\1\2||", s)
    s = re.sub(r"\|\\mathcal\{([A-Za-z0-9]+)\}_\{?([A-Za-z0-9]+)\}?\|", r"|\1\2|", s)
    s = re.sub(r"\\(?:mathrm|text|mathcal|mathbb|mathtt|texttt)\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\\(?:times|cdot|oplus|rtimes|mid)", "x", s)
    s = re.sub(r"\\(?:le|ge|neq|approx|sim|to|rightarrow|in)", "==", s)
    s = re.sub(r"\\(?:leftrightarrow)", "<->", s)
    s = re.sub(r"\\(?:dots|ldots|cdots)", "...", s)
    s = re.sub(r"\\(?:lceil|rceil|langle|rangle)", "|", s)
    s = re.sub(r"\\[A-Za-z]+", "M", s)
    s = s.replace("{", "").replace("}", "").replace("^", "").replace("_", "")
    return s


def cell_metrics(cell: str) -> tuple[float, float]:
    """Return (full_width_chars, max_unbreakable_chars) in approximate LaTeX width units."""
    # Strip markdown links [text](url) -> text
    s = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", cell)
    # Tokenize into backtick code `...`, math $...$, and plain text
    pattern = re.compile(r"(`[^`]+`|\$[^$]+\$)")
    parts = pattern.split(s)

    total_w = 0.0
    max_unbreak = 0.0
    cur_unbreak = 0.0

    def flush_space() -> None:
        nonlocal max_unbreak, cur_unbreak
        if cur_unbreak > max_unbreak:
            max_unbreak = cur_unbreak
        cur_unbreak = 0.0

    for part in parts:
        if not part:
            continue
        if part.startswith("`") and part.endswith("`"):
            code = part[1:-1]
            # Monospace is slightly wider than roman; breaks only at space or '_' (via \allowbreak)
            chunks = re.split(r"( +|_)", code)
            for ch in chunks:
                if not ch:
                    continue
                if ch.startswith(" "):
                    total_w += len(ch)
                    flush_space()
                elif ch == "_":
                    total_w += 1.12
                    cur_unbreak += 1.12
                    flush_space()
                else:
                    w = len(ch) * 1.12
                    total_w += w
                    cur_unbreak += w
        elif part.startswith("$") and part.endswith("$"):
            mt = approx_math_text(part[1:-1])
            # Math can break around binary/rel ops or spaces, approximate unbreakable chunks
            chunks = re.split(r"( +|=|\+|-|,)", mt)
            for ch in chunks:
                if not ch:
                    continue
                if ch.startswith(" "):
                    total_w += 0.8
                    flush_space()
                elif ch in ("=", "+", "-", ","):
                    total_w += 1.0
                    cur_unbreak += 1.0
                    flush_space()
                else:
                    w = len(ch) * 0.95
                    total_w += w
                    cur_unbreak += w
        else:
            plain = part.replace("**", "").replace("*", "")
            chunks = re.split(r"( +|-)", plain)
            for ch in chunks:
                if not ch:
                    continue
                if ch.startswith(" "):
                    total_w += len(ch) * 0.85
                    flush_space()
                elif ch == "-":
                    total_w += 0.9
                    cur_unbreak += 0.9
                    flush_space()
                else:
                    w = len(ch) * 1.00
                    total_w += w
                    cur_unbreak += w

    flush_space()
    return total_w, max_unbreak


def parse_align(tok: str) -> tuple[bool, bool, int]:
    t = tok.strip()
    left = t.startswith(":")
    right = t.endswith(":")
    min_len = (1 if left else 0) + (1 if right else 0) + 1
    min_len = max(min_len, 3)
    return left, right, min_len


def make_sep_token(left: bool, right: bool, length: int) -> str:
    dashes = length - (1 if left else 0) - (1 if right else 0)
    assert dashes >= 1, (left, right, length)
    return (":" if left else "") + ("-" * dashes) + (":" if right else "")


def allocate_tokens(rows: list[list[str]], aligns: list[tuple[bool, bool, int]], total_units: int = 80) -> list[int]:
    ncols = len(aligns)
    # Effective chars per unit out of 80:
    # Usable width = 576 - 12*(ncols - 1) pt; 1 monospace-equiv char ~ 5.0 pt => chars_per_unit = usable / (5.0 * 80)
    usable_pt = 576.0 - 12.0 * (ncols - 1)
    chars_per_unit = usable_pt / (5.05 * total_units)

    max_full = [0.0] * ncols
    max_unb = [0.0] * ncols
    for r in rows:
        for j in range(min(ncols, len(r))):
            fw, uw = cell_metrics(r[j])
            if fw > max_full[j]:
                max_full[j] = fw
            if uw > max_unb[j]:
                max_unb[j] = uw

    floors = []
    desires = []
    for j in range(ncols):
        _, _, min_tok = aligns[j]
        # Floor guarantees the longest unbreakable chunk fits in p{...} without Overfull \hbox
        unb_units = math.ceil(max_unb[j] / chars_per_unit + 0.65)
        fl = max(min_tok, unb_units)
        floors.append(fl)
        # Concave desire: full width up to 22 chars, then 0.50 slope above 22
        m = max_full[j]
        des_chars = min(m, 22.0) + 0.50 * max(0.0, m - 22.0)
        des_units = max(float(fl), des_chars / chars_per_unit)
        desires.append(des_units)

    sum_floors = sum(floors)
    if sum_floors >= total_units:
        # Scale down proportionally while keeping min_tok
        min_toks = [a[2] for a in aligns]
        rem = total_units - sum(min_toks)
        extra = [max(0, floors[j] - min_toks[j]) for j in range(ncols)]
        s_ex = sum(extra) or 1
        alloc = [min_toks[j] + int(math.floor(rem * extra[j] / s_ex)) for j in range(ncols)]
        while sum(alloc) < total_units:
            best = max(range(ncols), key=lambda j: (rem * extra[j] / s_ex) - (alloc[j] - min_toks[j]))
            alloc[best] += 1
        return alloc

    rem = total_units - sum_floors
    weights = [max(0.05, desires[j] - floors[j]) for j in range(ncols)]
    w_sum = sum(weights)
    raw_add = [rem * w / w_sum for w in weights]
    add_int = [int(math.floor(x)) for x in raw_add]
    alloc = [floors[j] + add_int[j] for j in range(ncols)]
    leftover = total_units - sum(alloc)
    fracs = sorted(range(ncols), key=lambda j: raw_add[j] - add_int[j], reverse=True)
    for j in fracs[:leftover]:
        alloc[j] += 1
    return alloc


def process_file(path: Path) -> int:
    lines = path.read_text().splitlines()
    in_code = False
    changed = 0
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("```"):
            in_code = not in_code
            i += 1
            continue
        if (
            not in_code
            and i + 1 < len(lines)
            and lines[i].lstrip().startswith("|")
            and SEP_RE.match(lines[i + 1].strip())
        ):
            start = i
            sep_idx = i + 1
            i += 2
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                i += 1
            end = i
            max_raw = max(len(lines[k]) for k in range(start, end))
            if max_raw <= 72:
                continue
            sep_cells = split_row(lines[sep_idx])
            aligns = [parse_align(c) for c in sep_cells]
            rows = [split_row(lines[start])] + [split_row(lines[k]) for k in range(sep_idx + 1, end)]
            alloc = allocate_tokens(rows, aligns, total_units=80)
            sep_tokens = [make_sep_token(aligns[j][0], aligns[j][1], alloc[j]) for j in range(len(aligns))]
            new_sep = "| " + " | ".join(sep_tokens) + " |"
            assert len(new_sep) <= 120, (path, start + 1, len(new_sep), new_sep)
            if lines[sep_idx] != new_sep:
                lines[sep_idx] = new_sep
                changed += 1
            continue
        i += 1
    if changed:
        path.write_text("\n".join(lines) + "\n")
    return changed


def main() -> int:
    targets = sys.argv[1:] or [str(ROOT / "benchmark" / "BENCHMARK.md"), str(ROOT / "benchmark" / "HISTORY.md")]
    for rel in targets:
        p = Path(rel)
        n = process_file(p)
        print(f"{rel}: updated {n} wide pipe tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
