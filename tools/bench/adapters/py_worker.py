#!/usr/bin/env python3
"""Fresh-process worker for the `dist_m4ri.py` Python adapter (`tools/bench/adapters/dist_m4ri_py.py`).

Imports `dist_m4ri` strictly from `spec['repo_root']`, disables the distance cache, calls the
appropriate `compute_*_distance(..., method=3, threads=..., timeout=..., seed=..., engine=...)`
API with default-league options only, and emits `@@PYRESULT {json}` on the final stdout line.
"""
from __future__ import annotations

import json
import os
import sys
import time
from typing import Any, Dict


def main() -> int:
    t0 = time.monotonic()
    spec = json.loads(sys.argv[1])
    # Remove this script's directory (tools/bench/adapters/) so adapters/qec_cert.py never shadows the qec_cert package.
    self_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path[:] = [p for p in sys.path if os.path.abspath(p or ".") != self_dir]
    repo_root = os.path.abspath(spec["repo_root"])
    sys.path.insert(0, repo_root)
    qec_root = spec.get("qec_cert_root")
    if qec_root and os.path.isdir(qec_root):
        qec_abs = os.path.abspath(qec_root)
        if qec_abs not in sys.path:
            sys.path.insert(1, qec_abs)

    t_imp0 = time.monotonic()
    import dist_m4ri as dm  # noqa: E402

    t_import = time.monotonic() - t_imp0
    dm_dir = os.path.dirname(os.path.abspath(dm.__file__))
    if dm_dir != repo_root:
        raise RuntimeError(f"dist_m4ri imported from {dm_dir}, expected {repo_root}")

    dm.disable_distance_cache()

    bin_path = spec.get("binary_path")
    if bin_path:
        bin_path = os.path.abspath(bin_path)

    engine = spec.get("engine")  # None (off), "auto", or "cert"
    kw: Dict[str, Any] = {
        "method": 3,
        "threads": int(spec["threads"]),
        "timeout": float(spec["budget_s"]),
        "seed": int(spec["seed"]),
        "debug": int(spec.get("debug", 0)),
    }
    if bin_path:
        kw["dist_m4ri"] = bin_path
    if engine is not None:
        kw["engine"] = engine

    out: Dict[str, Any] = {
        "dm_version": getattr(dm, "__version__", None),
        "bin_path": bin_path or dm.find_dist_m4ri_binary(),
        "engine": engine,
        "budget_s": float(spec["budget_s"]),
        "dist": None,
        "d_info": None,
        "error": None,
    }

    t_call0 = time.monotonic()
    try:
        kind = spec["kind"]
        files = spec["files"]
        if kind == "css_sector":
            res = dm.compute_quantum_distance(files["H"], G=files["G"], return_info=True, **kw)
            out["dist"], out["d_info"] = int(res[0]), [int(x) for x in res[1]]
        elif kind == "hl":
            res = dm.compute_quantum_distance(files["H"], L=files["L"], return_info=True, **kw)
            out["dist"], out["d_info"] = int(res[0]), [int(x) for x in res[1]]
        elif kind == "classical":
            res = dm.compute_classical_distance(files["H"], return_info=True, **kw)
            out["dist"], out["d_info"] = int(res[0]), [int(x) for x in res[1]]
        elif kind == "dem":
            res = dm.compute_dem_distance(dem=files["dem"], **kw)
            out["dist"], out["d_info"] = int(res[0]), [int(x) for x in res[1]]
        elif kind == "circuit":
            res = dm.compute_dem_distance(circuit=files["circuit"], **kw)
            out["dist"], out["d_info"] = int(res[0]), [int(x) for x in res[1]]
        else:
            raise ValueError(f"Unsupported kind {kind!r}")
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"[:400]
    t_call = time.monotonic() - t_call0

    stats = getattr(dm, "_last_run_stats", None) or {}
    out.update(
        t_import=round(t_import, 6),
        t_call=round(t_call, 6),
        t_proc=round(time.monotonic() - t0, 6),
        qec_imported="qec_cert" in sys.modules,
        last_run_stats=json.loads(json.dumps(stats, default=str)),
    )
    print("@@PYRESULT " + json.dumps(out), flush=True)
    return 0 if out["error"] is None else 1


if __name__ == "__main__":
    raise SystemExit(main())
