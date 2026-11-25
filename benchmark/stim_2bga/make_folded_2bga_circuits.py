#!/usr/bin/env python3
"""Deterministic generator for Shaw-Terhal folded sequential 2BGA Stim circuits.

Reproduces the noiseless folded memory circuit templates (T = 1 and T = d, X and Z bases) from L. P. Pryadko's
hr-2bga optcirc campaign (following Shaw & Terhal, arXiv:2604.09797):
  - 2bga_162_24_16_C9xC9_folded_T{1,16}_{X,Z}.stim: [[162, 24, 16]] over C9 x C9 (w = 10 = 5 + 5, d_circ in [8, 9])
  - 2bga_254_44_19_C127_folded_T{1,19}_{X,Z}.stim:  [[254, 44, 19]] over C127 (w = 12 = 6 + 6, d_circ in [9, 17])

Usage:
  python3 benchmark/stim_2bga/make_folded_2bga_circuits.py --out-dir scratch/folded_2bga [--only 2bga_254]
"""
from __future__ import annotations

import argparse
from pathlib import Path
import stim

# GAP AsSSortedList(SmallGroup(81, 2)) coordinates (i, j) in Z_9 x Z_9 for g1^i * g2^j:
PAIRS_81_2 = [
    (0, 0), (1, 0), (0, 1), (3, 0), (0, 3), (2, 0), (1, 1), (4, 0), (1, 3),
    (0, 2), (3, 1), (0, 4), (6, 0), (3, 3), (0, 6), (2, 1), (5, 0), (2, 3),
    (1, 2), (4, 1), (1, 4), (7, 0), (4, 3), (1, 6), (3, 2), (0, 5), (6, 1),
    (3, 4), (0, 7), (6, 3), (3, 6), (2, 2), (5, 1), (2, 4), (8, 0), (5, 3),
    (2, 6), (4, 2), (1, 5), (7, 1), (4, 4), (1, 7), (7, 3), (4, 6), (6, 2),
    (3, 5), (0, 8), (6, 4), (3, 7), (6, 6), (5, 2), (2, 5), (8, 1), (5, 4),
    (2, 7), (8, 3), (5, 6), (7, 2), (4, 5), (1, 8), (7, 4), (4, 7), (7, 6),
    (6, 5), (3, 8), (6, 7), (8, 2), (5, 5), (2, 8), (8, 4), (5, 7), (8, 6),
    (7, 5), (4, 8), (7, 7), (6, 8), (8, 5), (5, 8), (8, 7), (7, 8), (8, 8),
]


def bits(x: int) -> list[int]:
    out = []
    while x:
        low = x & -x
        out.append(low.bit_length() - 1)
        x ^= low
    return out


def echelon(rows: list[int]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for r in rows:
        while r:
            p = r.bit_length() - 1
            q = basis.get(p)
            if q is None:
                basis[p] = r
                break
            r ^= q
    return basis


def logical_basis(h_rows: list[int], g_rows: list[int], n: int) -> list[int]:
    b_h = echelon(h_rows)
    piv = sorted(b_h)
    for idx, p in enumerate(piv):
        r = b_h[p]
        for q in piv[idx + 1:]:
            if (b_h[q] >> p) & 1:
                b_h[q] ^= r
    ker = []
    for f in range(n):
        if f in b_h:
            continue
        x = 1 << f
        for p, r in b_h.items():
            if (r >> f) & 1:
                x |= 1 << p
        ker.append(x)
    b_g = echelon(g_rows)
    out = []
    for x in ker:
        r = x
        while r:
            p = r.bit_length() - 1
            q = b_g.get(p)
            if q is None:
                b_g[p] = r
                out.append(x)
                break
            r ^= q
    return out


def group_tables_c9xc9() -> tuple[list[list[int]], list[int], list[tuple[int, int]]]:
    idx_of = {pt: k for k, pt in enumerate(PAIRS_81_2)}
    mul = [[idx_of[((i1 + j1) % 9, (i2 + j2) % 9)] for (j1, j2) in PAIRS_81_2] for (i1, i2) in PAIRS_81_2]
    inv = [idx_of[((-i1) % 9, (-i2) % 9)] for (i1, i2) in PAIRS_81_2]
    return mul, inv, PAIRS_81_2


def group_tables_cyclic(ell: int) -> tuple[list[list[int]], list[int], list[tuple[int, int]]]:
    mul = [[(i + j) % ell for j in range(ell)] for i in range(ell)]
    inv = [(-i) % ell for i in range(ell)]
    coords_ij = [(i, 0) for i in range(ell)]
    return mul, inv, coords_ij


def build_folded_2bga(
    mul: list[list[int]],
    inv: list[int],
    coords_ij: list[tuple[int, int]],
    a_1based: list[int],
    b_1based: list[int],
    order_terms: list[int],
    basis: str,
    rounds: int,
) -> stim.Circuit:
    ell = len(mul)
    a = [x - 1 for x in a_1based]
    b = [x - 1 for x in b_1based]
    wa, wb = len(a), len(b)
    w = wa + wb
    xq = [[mul[inv[g]][al] for g in a] + [ell + mul[al][inv[h]] for h in b] for al in range(ell)]
    zq = [[ell + mul[g][be] for g in a] + [mul[be][h] for h in b] for be in range(ell)]
    hx = [sum(1 << q for q in row) for row in xq]
    hz = [sum(1 << q for q in row) for row in zq]

    xa = [2 * ell + al for al in range(ell)]
    za = [3 * ell + be for be in range(ell)]
    data = list(range(2 * ell))

    coords: list[tuple[int, int, int]] = [(0, 0, 0)] * (4 * ell)
    for al, (i, j) in enumerate(coords_ij):
        coords[al] = (2 * i + 1, 2 * j, 0)
        coords[ell + al] = (2 * i, 2 * j + 1, 0)
        coords[xa[al]] = (2 * i, 2 * j, 1)
        coords[za[al]] = (2 * i + 1, 2 * j + 1, 2)

    def layer_x(ti: int) -> list[int]:
        return [v for al in range(ell) for v in (xa[al], xq[al][ti])]

    def layer_z(ti: int) -> list[int]:
        return [v for be in range(ell) for v in (zq[be][ti], za[be])]

    def round_body(backward: bool) -> tuple[stim.Circuit, list[tuple[str, int]]]:
        body = stim.Circuit()
        blocks = [("X", order_terms), ("Z", order_terms)]
        meas_all = [("X", al) for al in range(ell)] + [("Z", be) for be in range(ell)]
        if backward:
            for typ, ord_t in blocks[::-1]:
                body.append("RX" if typ == "X" else "R", xa if typ == "X" else za)
                body.append("TICK")
                for ti in ord_t[::-1]:
                    body.append("CX", layer_x(ti) if typ == "X" else layer_z(ti))
                    body.append("TICK")
            body.append("MX", xa)
            body.append("M", za)
            body.append("TICK")
            return body, meas_all
        body.append("RX", xa)
        body.append("R", za)
        body.append("TICK")
        meas: list[tuple[str, int]] = []
        for typ, ord_t in blocks:
            for ti in ord_t:
                body.append("CX", layer_x(ti) if typ == "X" else layer_z(ti))
                body.append("TICK")
            body.append("MX" if typ == "X" else "M", xa if typ == "X" else za)
            body.append("TICK")
            meas += [(typ, al) for al in range(ell)]
        return body, meas

    bodies = {"F": round_body(False), "B": round_body(True)}
    m_cnt = len(bodies["F"][1])
    anc_xy = {("X", al): coords[xa[al]][:2] for al in range(ell)}
    anc_xy.update({("Z", be): coords[za[be]][:2] for be in range(ell)})

    def detectors(cur: str, prev: str | None) -> stim.Circuit:
        d = stim.Circuit()
        pos = {m: idx for idx, m in enumerate(bodies[cur][1])}
        ppos = {m: idx for idx, m in enumerate(bodies[prev][1])} if prev else {}
        for m in bodies[cur][1]:
            if prev is None and m[0] != basis:
                continue
            tg = [stim.target_rec(pos[m] - m_cnt)] + ([stim.target_rec(ppos[m] - 2 * m_cnt)] if prev else [])
            d.append("DETECTOR", tg, [anc_xy[m][0], anc_xy[m][1], 0])
        d.append("SHIFT_COORDS", [], [0, 0, 1])
        return d

    def segment(cur: str, prev: str | None) -> stim.Circuit:
        return bodies[cur][0] + detectors(cur, prev)

    plan = ["B" if r % 2 == 1 else "F" for r in range(rounds)]
    c = stim.Circuit()
    for q, (x, y, typ) in enumerate(coords):
        c.append("QUBIT_COORDS", [q], [x, y, typ])
    c.append("R" if basis == "Z" else "RX", data)
    c.append("TICK")
    c += segment("F", None)
    pairs = (len(plan) - 1) // 2
    if pairs:
        c.append(stim.CircuitRepeatBlock(pairs, segment("B", "F") + segment("F", "B")))
    if len(plan) % 2 == 0:
        c += segment("B", "F")
    c.append("M" if basis == "Z" else "MX", data)
    nd = len(data)
    last = {m: idx for idx, m in enumerate(bodies[plan[-1]][1])}
    for idx in range(ell):
        m = (basis, idx)
        supp = zq[idx] if basis == "Z" else xq[idx]
        tg = [stim.target_rec(last[m] - m_cnt - nd)] + [stim.target_rec(q - nd) for q in sorted(supp)]
        c.append("DETECTOR", tg, [anc_xy[m][0], anc_xy[m][1], 0])
    logs = logical_basis(hx, hz, 2 * ell) if basis == "Z" else logical_basis(hz, hx, 2 * ell)
    assert len(logs) == 2 * ell - len(echelon(hx)) - len(echelon(hz))
    for idx, x in enumerate(logs):
        c.append("OBSERVABLE_INCLUDE", [stim.target_rec(q - nd) for q in bits(x)], idx)
    return c


SPECS = {
    "2bga_162": {
        "prefix": "2bga_162_24_16_C9xC9_folded",
        "group": group_tables_c9xc9,
        "a": [1, 17, 31, 37, 42],
        "b": [1, 7, 39, 59, 72],
        "order": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        "rounds": [1, 16],
    },
    "2bga_254": {
        "prefix": "2bga_254_44_19_C127_folded",
        "group": lambda: group_tables_cyclic(127),
        "a": [1, 18, 25, 39, 98, 125],
        "b": [1, 26, 36, 44, 72, 76],
        # a4,a5,a3,a2,a1,a0,b0,b1,b2,b4,b3,b5:
        "order": [4, 5, 3, 2, 1, 0, 6, 7, 8, 10, 9, 11],
        "rounds": [1, 19],
    },
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True, help="Output directory for generated .stim files")
    ap.add_argument("--only", choices=sorted(SPECS), help="Generate only one code family (default: both)")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    targets = [args.only] if args.only else list(SPECS)
    for key in targets:
        sp = SPECS[key]
        mul, inv, cij = sp["group"]()
        for t_val in sp["rounds"]:
            for bas in ("X", "Z"):
                circ = build_folded_2bga(mul, inv, cij, sp["a"], sp["b"], sp["order"], bas, t_val)
                out_p = args.out_dir / f"{sp['prefix']}_T{t_val}_{bas}.stim"
                out_p.write_text(str(circ) + "\n", encoding="utf-8")
                print(
                    f"Wrote {out_p} (qubits={circ.num_qubits}, dets={circ.num_detectors}, obs={circ.num_observables})"
                )


if __name__ == "__main__":
    main()
