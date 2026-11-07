#!/usr/bin/env python3
"""Memory-experiment Stim circuits for the mitten codes of a7b/yarn (arXiv:2607.28795), with the published hook-free
syndrome-extraction schedule (processor_codes/<family>/[[n,k,d]]/hook_free_SE_cycle_schedule.json).

Each cycle: reset the X- and Z-check ancillas; the 12 layers of the X-check round (CX from the X-check ancilla to the
data qubit for every listed pair [check row of Hx, data column]); MX of the X-check ancillas; the 12 layers of the
Z-check round (CX from the data qubit to the Z-check ancilla); M of the Z-check ancillas. The cycles after the first are
in a REPEAT block. The data qubits start in |0> (Z basis) or |+> (X basis) and are measured in that basis at the end.
Detectors: the checks of the memory basis in the first cycle, every check in the later cycles (compared with the
previous cycle), and the checks of the memory basis from the final data measurements; observables: the rows of Lz
(Z basis) or Lx (X basis). Noise: uniform circuit-level noise of dist_m4ri.add_noise(p) (DEPOLARIZE2 after each CX,
DEPOLARIZE1(p/10) on idle qubits at each TICK, flips of resets and measurements).

Usage, from a checkout of https://github.com/a7b/yarn (commit 1b58edd24e19b6b5da6814fd3a74ebb376b6e90f):
    python3 benchmark/stim_mitten/make_yarn_mitten_circuits.py "YARN/processor_codes/mitten/[[150,30,10]]" \
        benchmark/stim_mitten [--rounds R] [--p P] [--basis X|Z|XZ] [--check-rounds R]
(R = 0: the d of the code name). Output: yarn_<family>_<n>_<k>_<d>_hookfree_r<R>_<basis>.stim.
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import stim

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
import dist_m4ri  # noqa: E402  (add_noise)


def load(code_dir):
    m = {k: np.load(os.path.join(code_dir, k + ".npy")).astype(np.uint8) for k in ("Hx", "Hz", "Lx", "Lz")}
    sched = json.load(open(os.path.join(code_dir, "hook_free_SE_cycle_schedule.json")))
    return m, sched


def layers_of(sched, key, h, nq_data):
    """The layers of one check round as lists of (check row, data column); checks against h and qubit conflicts."""
    layers, seen = [], np.zeros_like(h)
    for layer in sched[key]["layers"]:
        pairs = [(int(i), int(j)) for mv in layer["moves"] for i, j in mv["cz_gates"]]
        checks = [i for i, _ in pairs]
        datas = [j for _, j in pairs]
        assert len(set(checks)) == len(checks) and len(set(datas)) == len(datas), f"qubit used twice in a layer"
        for i, j in pairs:
            assert h[i, j] == 1, f"{key}: pair ({i}, {j}) is not in the check matrix"
            assert seen[i, j] == 0, f"{key}: pair ({i}, {j}) listed twice"
            seen[i, j] = 1
        layers.append(pairs)
    assert (seen == h).all(), f"{key}: the layers do not cover the check matrix"
    return layers


def build(m, sched, basis, rounds):
    """The noiseless memory circuit: the first cycle, then the other rounds - 1 cycles in a REPEAT block."""
    hx, hz = m["Hx"], m["Hz"]
    n, mx, mz = hx.shape[1], hx.shape[0], hz.shape[0]
    xl = layers_of(sched, "X_layers", hx, n)
    zl = layers_of(sched, "Z_layers", hz, n)
    xa = [n + i for i in range(mx)]
    za = [n + mx + i for i in range(mz)]
    mt = mx + mz  # measurements per cycle: X checks (lookbacks -mt .. -mz-1), then Z checks (-mz .. -1)

    def cycle(first):
        c = stim.Circuit()
        c.append("RX", xa)
        c.append("R", za)
        c.append("TICK")
        for pairs in xl:
            c.append("CX", [q for i, j in pairs for q in (xa[i], j)])
            c.append("TICK")
        c.append("MX", xa)
        c.append("TICK")
        for pairs in zl:
            c.append("CX", [q for i, j in pairs for q in (j, za[i])])
            c.append("TICK")
        c.append("M", za)
        for typ, num, off in (("X", mx, -mt), ("Z", mz, -mz)):
            if first and typ != basis:
                continue
            for i in range(num):
                targets = [stim.target_rec(off + i)] + ([] if first else [stim.target_rec(off + i - mt)])
                c.append("DETECTOR", targets, [i, 0, 0 if typ == "X" else 1])
        c.append("SHIFT_COORDS", [], [0, 1])
        c.append("TICK")
        return c

    c = stim.Circuit()
    c.append("RX" if basis == "X" else "R", range(n))
    c += cycle(True)
    if rounds > 1:
        c.append(stim.CircuitRepeatBlock(rounds - 1, cycle(False)))
    c.append("MX" if basis == "X" else "M", range(n))
    h, lo = (hx, m["Lx"]) if basis == "X" else (hz, m["Lz"])
    last = -n - mt if basis == "X" else -n - mz  # lookback of check 0 of the memory basis in the last cycle
    for i in range(h.shape[0]):
        targets = [stim.target_rec(-n + j) for j in np.nonzero(h[i])[0]] + [stim.target_rec(last + i)]
        c.append("DETECTOR", targets, [i, 0, 0 if basis == "X" else 1])
    for k in range(lo.shape[0]):
        c.append("OBSERVABLE_INCLUDE", [stim.target_rec(-n + j) for j in np.nonzero(lo[k])[0]], k)
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("code_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--rounds", type=int, default=0)
    ap.add_argument("--p", type=float, default=0.001)
    ap.add_argument("--basis", default="XZ")
    ap.add_argument("--check-rounds", type=int, default=2, help="rounds of the noisy test circuit for the DEM check")
    a = ap.parse_args()
    m, sched = load(a.code_dir)
    nkd = [int(x) for x in re.findall(r"\d+", os.path.basename(os.path.normpath(a.code_dir)))]
    fam = os.path.basename(os.path.dirname(os.path.normpath(a.code_dir)))
    prefix = {"mitten": "mitten", "structured_mitten": "smitten"}.get(fam, fam)
    rounds = a.rounds or nkd[2]
    os.makedirs(a.out_dir, exist_ok=True)
    for basis in a.basis:
        # determinism of all detectors and observables: stim raises otherwise
        test = dist_m4ri.add_noise(build(m, sched, basis, a.check_rounds), p=a.p)
        dem = test.detector_error_model(decompose_errors=False)
        c = dist_m4ri.add_noise(build(m, sched, basis, rounds), p=a.p)
        name = f"yarn_{prefix}_{nkd[0]}_{nkd[1]}_{nkd[2]}_hookfree_r{rounds}_{basis}.stim"
        header = (f"# [[{nkd[0]},{nkd[1]},{nkd[2]}]] {fam} code of a7b/yarn (arXiv:2607.28795), group "
                  f"{sched['code'].get('group')}; {basis}-basis memory, {rounds} rounds, published hook-free SE "
                  f"schedule (X-check round, then Z-check round, 12 CX layers each); dist_m4ri.add_noise(p={a.p})\n")
        with open(os.path.join(a.out_dir, name), "w") as f:
            f.write(header + str(c) + "\n")
        print(json.dumps({"file": name, "qubits": c.num_qubits, "detectors": c.num_detectors,
                          "observables": c.num_observables, f"dem_errors_r{a.check_rounds}": dem.num_errors}),
              flush=True)


if __name__ == "__main__":
    main()
