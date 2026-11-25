#!/usr/bin/env python3
"""Deterministic generator for large (> 2 MB) A != 0 detector-drop DEM benchmarks.

Generates <prefix>_dropA.dem and <prefix>_dropA.drop from benchmark/stim_bb/ circuits by perturbing
logical observables as L(j) ^= A * H_M(j) on the minority-basis detector set M with fixed seed=12345:
  - dgross_288_dropA.{dem,drop}: from double_gross_288_12_18_si1000_r18_X.stim (~17 MB DEM, |M|=2448, k=12)
  - bb_756_dropA.{dem,drop}:     from bb_756_16_34_si1000_r34_Z.stim (> 50 MB DEM, |M|=12474, k=16)

Usage:
  python3 benchmark/dem/make_large_drop_dems.py --out-dir scratch/large_drop [--only dgross_288]
"""
from __future__ import annotations

import argparse
import random
import subprocess
from pathlib import Path
import stim

SPECS = {
    "dgross_288": "double_gross_288_12_18_si1000_r18_X.stim",
    "bb_756": "bb_756_16_34_si1000_r34_Z.stim",
}


def generate_one(repo_root: Path, name: str, stim_file: str, out_dir: Path, seed: int = 12345) -> None:
    bin_path = repo_root / "src" / "dist_m4ri"
    stim_path = repo_root / "benchmark" / "stim_bb" / stim_file
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_src = out_dir / f"_{name}_src.dem"
    drop_path = out_dir / f"{name}_dropA.drop"
    out_dem = out_dir / f"{name}_dropA.dem"

    circ = stim.Circuit.from_file(str(stim_path))
    dem = circ.detector_error_model(decompose_errors=False, flatten_loops=True)
    dem.to_file(str(tmp_src))

    subprocess.run(
        [str(bin_path), f"fdem={tmp_src}", f"outD={drop_path}", "method=2", "wmax=1", "threads=4", "debug=1"],
        check=True, capture_output=True, text=True,
    )
    dropped: list[int] = []
    for ln in drop_path.read_text(encoding="utf-8").splitlines():
        if ln.strip() and not ln.startswith("#"):
            dropped.extend(int(x) for x in ln.split())
    drop_idx = {d: q for q, d in enumerate(dropped)}
    num_m = len(dropped)
    k = dem.num_observables
    rng = random.Random(seed)
    a_mask = [rng.randrange(1, 1 << k) if q == 0 else rng.randrange(0, 1 << k) for q in range(num_m)]

    out_lines = [
        f"# Perturbed DEM from {stim_file}: L(j) ^= A * H_M(j) (|M|={num_m}, k={k}, seed={seed})",
    ]
    for ln in tmp_src.read_text(encoding="utf-8").splitlines():
        if not ln.startswith("error"):
            out_lines.append(ln)
            continue
        toks = ln.split()
        head = toks[0]
        dets = []
        lmask = 0
        hm_xor = 0
        for tok in toks[1:]:
            if tok == "^":
                continue
            if tok.startswith("D"):
                d = int(tok[1:])
                dets.append(tok)
                q = drop_idx.get(d)
                if q is not None:
                    hm_xor ^= a_mask[q]
            elif tok.startswith("L"):
                lmask ^= 1 << int(tok[1:])
        new_lmask = lmask ^ hm_xor
        ltoks = [f"L{i}" for i in range(k) if (new_lmask >> i) & 1]
        out_lines.append(" ".join([head] + dets + ltoks))

    out_dem.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    tmp_src.unlink()
    print(f"Wrote {out_dem} and {drop_path} (|M|={num_m}, k={k}, seed={seed})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True, help="Output directory for generated .dem and .drop files")
    ap.add_argument("--only", choices=sorted(SPECS), help="Generate only one target (default: both)")
    ap.add_argument("--seed", type=int, default=12345, help="Fixed RNG seed for A != 0 (default: 12345)")
    args = ap.parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    targets = [args.only] if args.only else list(SPECS)
    for name in targets:
        generate_one(repo_root, name, SPECS[name], args.out_dir, seed=args.seed)


if __name__ == "__main__":
    main()
