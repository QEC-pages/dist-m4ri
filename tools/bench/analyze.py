#!/usr/bin/env python3
"""
tools/bench/analyze.py: Pure analysis, tables, figures, and report generator for the
dist-m4ri certification benchmark (tools/bench/SCHEMA.md, dm4ri-bench-run/1).

Public API:
- `env_id(env) -> str`: canonical sha1 digest of an env.json dict.
- `load_session(session_dir) -> tuple[list[dict], dict | None, dict]`: load one session directory.
- `load_sessions(session_dirs) -> tuple[list[dict], dict[str, dict], list[dict]]`: load multiple sessions.
- `analyze(manifest, records, envs=None, *, baseline=..., exclude_contaminated=False, n_boot=2000,
  boot_seed=1) -> dict`: pure, deterministic, JSON-serializable analysis (Tables A-D, Figures 1-5 +
  optional canary drift, and all captions). No filesystem I/O.
- `render_figures(result, out_dir) -> dict[str, str]`: render PNG figures with matplotlib (Agg) from
  `result["figures"]` data only.
- `render_report(result, fig_files=None, *, load_meta=None) -> str`: format REPORT.md (all lines <= 120
  columns, zero hand-written verdicts).

CLI usage:
    python3 tools/bench/analyze.py --manifest benchmark/manifest/certify_v1.json \\
        --session path/to/session1 [--session path/to/session2 ...] \\
        --out path/to/out [--exclude-contaminated] [--baseline dist_m4ri:default@p128] [--no-figures]
"""
from __future__ import annotations

import argparse
import atexit
import hashlib
import json
import math
import os
import shutil
import tempfile
import textwrap
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

SUMMARY_SCHEMA = "dm4ri-bench-summary/1"
DEFAULT_BASELINE = "dist_m4ri:default@p128"
TIER_ORDER_DEFAULT = ("T0", "T1", "T2", "T3")
VALID_BRACKET_STATUSES = frozenset({"certified", "bracket", "timeout"})
WALL_TIME_SENTENCE = (
    "Wall time is process wall time from fork to exit (measured by the harness), "
    "including start-up, input stage, and Python imports."
)

_MPL_TMPDIR: str | None = None
_PLT_MOD: Any = None


# ======================================================================================================================
# Helpers: canonical hashing, formatting, JSON sanitization, deterministic RNG
# ======================================================================================================================

def env_id(env: Mapping[str, Any]) -> str:
    """Return the sha1 hex digest of the canonical JSON representation of an env.json mapping."""
    payload = json.dumps(env, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()


def _entry_name(rec: Mapping[str, Any]) -> str:
    """Build the canonical engine entry identifier `name:config@p<threads>` for a run record."""
    eng = rec.get("engine") or {}
    opt = rec.get("options") or {}
    name = str(eng.get("name") or "unknown")
    cfg = str(eng.get("config") or "default")
    threads = opt.get("threads")
    t_str = str(int(threads)) if threads is not None else "?"
    return f"{name}:{cfg}@p{t_str}"


def _is_solved(rec: Mapping[str, Any]) -> bool:
    """Return True iff the record's final status is `certified` with a finite wall_sec."""
    res = rec.get("result") or {}
    if res.get("status") != "certified":
        return False
    w = _as_float(res.get("wall_sec"))
    return w is not None and w >= 0.0


def _as_float(val: Any) -> float | None:
    """Convert a value to a finite Python float, or None if missing/non-finite."""
    if val is None or isinstance(val, bool):
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _as_int(val: Any) -> int | None:
    """Convert a value to a Python int, or None if missing/invalid."""
    if val is None or isinstance(val, bool):
        return None
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def _record_budget(rec: Mapping[str, Any], tier_budgets: Mapping[str, float], tier: str) -> float:
    """Resolve the time budget T (seconds) for a run record."""
    opt = rec.get("options") or {}
    res = rec.get("result") or {}
    for cand in (opt.get("budget_s"), res.get("cap_sec"), tier_budgets.get(tier)):
        f = _as_float(cand)
        if f is not None and f > 0.0:
            return f
    return 0.0


def _record_par2(rec: Mapping[str, Any], tier_budgets: Mapping[str, float], tier: str) -> float:
    """Compute PAR-2 score: wall_sec if certified, else 2 * budget_s."""
    if _is_solved(rec):
        w = _as_float((rec.get("result") or {}).get("wall_sec"))
        if w is not None:
            return w
    return 2.0 * _record_budget(rec, tier_budgets, tier)


def _record_par10(rec: Mapping[str, Any], tier_budgets: Mapping[str, float], tier: str) -> float:
    """Compute PAR-10 score: wall_sec if certified, else 10 * budget_s."""
    if _is_solved(rec):
        w = _as_float((rec.get("result") or {}).get("wall_sec"))
        if w is not None:
            return w
    return 10.0 * _record_budget(rec, tier_budgets, tier)


def _stable_seed(base_seed: int, *parts: str) -> int:
    """Derive a deterministic 64-bit RNG seed from base_seed and label parts."""
    h = hashlib.sha256(("|".join([str(base_seed), *parts])).encode("utf-8")).digest()
    return int.from_bytes(h[:8], "little")


def _sanitize_json(obj: Any) -> Any:
    """Recursively convert numpy scalars and NaN/Inf floats into strict JSON-safe Python values."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return [_sanitize_json(x) for x in obj.tolist()]
    if isinstance(obj, Mapping):
        return {str(k): _sanitize_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_json(x) for x in obj]
    return str(obj)


def _fmt_time(val: float | None) -> str:
    """Format a duration in seconds compactly for Markdown tables."""
    if val is None or not math.isfinite(val):
        return "-"
    if val == 0.0:
        return "0"
    if val < 0.01:
        return f"{val:.4f}"
    if val < 10.0:
        return f"{val:.2f}"
    if val < 100.0:
        return f"{val:.1f}"
    return f"{val:.0f}"


def _fmt_ratio(val: float | None) -> str:
    """Format a ratio or speedup compactly."""
    if val is None or not math.isfinite(val):
        return "-"
    if val >= 100.0:
        return f"{val:.0f}x"
    if val >= 10.0:
        return f"{val:.1f}x"
    return f"{val:.2f}x"


def _fmt_num(val: float | None, decimals: int = 1) -> str:
    """Format a plain floating-point number or '-' if None."""
    if val is None or not math.isfinite(val):
        return "-"
    return f"{val:.{decimals}f}"


def _fmt_ref(d_ref: int | None, bracket_ref: Sequence[int] | None) -> str:
    """Format manifest reference distance or bracket."""
    if d_ref is not None:
        return str(int(d_ref))
    if bracket_ref and len(bracket_ref) == 2:
        lo = "-" if bracket_ref[0] is None else str(int(bracket_ref[0]))
        hi = "-" if bracket_ref[1] is None else str(int(bracket_ref[1]))
        return f"[{lo}, {hi}]"
    return "-"


def _short_inst_label(inst_id: str, max_len: int = 46) -> str:
    """Return an instance label trimmed cleanly if longer than max_len."""
    if len(inst_id) <= max_len:
        return inst_id
    return inst_id[: max_len - 3] + "..."


def _sort_entries(entries: Sequence[str], baseline: str) -> list[str]:
    """Sort engine entries with baseline first, then lexicographically."""
    unique = sorted(set(entries))
    if baseline in unique:
        unique.remove(baseline)
        return [baseline] + unique
    return unique


# ======================================================================================================================
# Session loading (thin I/O boundary)
# ======================================================================================================================

def load_session(session_dir: str | Path) -> tuple[list[dict[str, Any]], dict[str, Any] | None, dict[str, Any]]:
    """
    Load one benchmark session directory.

    Reads `results.jsonl` (skipping and counting malformed JSON lines) or falls back to
    `runs/*/run.json` when `results.jsonl` is absent. Optionally loads `env.json`, `plan.json`,
    and `session.json`.

    Args:
        session_dir: Path to a session directory written by `harness.py`.

    Returns:
        Tuple `(records, env, meta)` where `records` is the list of parsed run dicts in file order,
        `env` is the parsed `env.json` dict (or None), and `meta` contains loader metadata.
    """
    sdir = Path(session_dir)
    records: list[dict[str, Any]] = []
    bad_lines = 0
    bad_details: list[str] = []
    source = "empty"

    jsonl_path = sdir / "results.jsonl"
    if jsonl_path.is_file():
        source = "results.jsonl"
        with jsonl_path.open("r", encoding="utf-8", errors="replace") as fh:
            for lineno, raw_line in enumerate(fh, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except ValueError as exc:
                    bad_lines += 1
                    bad_details.append(f"{jsonl_path.name}:{lineno}: {exc}")
                    continue
                if isinstance(obj, dict):
                    records.append(obj)
                else:
                    bad_lines += 1
                    bad_details.append(f"{jsonl_path.name}:{lineno}: JSON value is not an object")
    else:
        runs_dir = sdir / "runs"
        if runs_dir.is_dir():
            source = "runs"
            for rpath in sorted(runs_dir.glob("*/run.json")):
                try:
                    obj = json.loads(rpath.read_text(encoding="utf-8", errors="replace"))
                except (OSError, ValueError) as exc:
                    bad_lines += 1
                    bad_details.append(f"{rpath.parent.name}/run.json: {exc}")
                    continue
                if isinstance(obj, dict):
                    records.append(obj)
                else:
                    bad_lines += 1
                    bad_details.append(f"{rpath.parent.name}/run.json: JSON value is not an object")

    def _opt_json(filename: str) -> dict[str, Any] | None:
        p = sdir / filename
        if not p.is_file():
            return None
        try:
            data = json.loads(p.read_text(encoding="utf-8", errors="replace"))
            return data if isinstance(data, dict) else None
        except (OSError, ValueError):
            bad_details.append(f"{filename}: unreadable JSON")
            return None

    env = _opt_json("env.json")
    plan = _opt_json("plan.json")
    session_info = _opt_json("session.json")
    meta = {
        "dir": str(sdir),
        "source": source,
        "records_loaded": len(records),
        "bad_lines": bad_lines,
        "bad_details": bad_details,
        "has_env": env is not None,
        "env_id": env_id(env) if env is not None else None,
        "has_plan": plan is not None,
        "has_session": session_info is not None,
        "session_id": (session_info or {}).get("session") or sdir.name,
    }
    return records, env, meta


def load_sessions(
    session_dirs: Sequence[str | Path],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Load multiple session directories in order, concatenating records and indexing envs by env_id."""
    all_records: list[dict[str, Any]] = []
    envs: dict[str, dict[str, Any]] = {}
    metas: list[dict[str, Any]] = []
    for sdir in session_dirs:
        recs, env, meta = load_session(sdir)
        all_records.extend(recs)
        if env is not None:
            eid = env_id(env)
            envs[eid] = env
        envs_dir = Path(sdir) / "envs"
        if envs_dir.is_dir():
            for ep in sorted(envs_dir.glob("*.json")):
                try:
                    eobj = json.loads(ep.read_text(encoding="utf-8"))
                    if isinstance(eobj, dict):
                        envs[env_id(eobj)] = eobj
                except (OSError, json.JSONDecodeError):
                    pass
        metas.append(meta)
    return all_records, envs, metas


# ======================================================================================================================
# Pure core analysis
# ======================================================================================================================

def _normalize_envs(
    envs: Mapping[str, dict[str, Any]] | Sequence[dict[str, Any]] | None,
) -> dict[str, dict[str, Any]]:
    """Normalize `envs` (mapping, sequence of env dicts, or None) to `{env_id: env_dict}`."""
    if envs is None:
        return {}
    if isinstance(envs, Mapping):
        out: dict[str, dict[str, Any]] = {}
        for k, v in envs.items():
            if isinstance(v, Mapping):
                key = str(k) if k else env_id(v)
                out[key] = dict(v)
        return out
    out = {}
    for item in envs:
        if isinstance(item, Mapping):
            out[env_id(item)] = dict(item)
    return out


def _build_manifest_index(
    manifest: Mapping[str, Any],
) -> tuple[list[str], dict[str, float], dict[str, dict[str, Any]], dict[str, list[str]]]:
    """Extract ordered tiers, tier budgets, per-instance metadata, and manifest-ordered tier lists."""
    tiers_raw = manifest.get("tiers") or {}
    tiers_order: list[str] = [t for t in TIER_ORDER_DEFAULT if t in tiers_raw]
    for t in tiers_raw:
        if t not in tiers_order:
            tiers_order.append(str(t))
    if not tiers_order:
        tiers_order = list(TIER_ORDER_DEFAULT)

    tier_budgets: dict[str, float] = {}
    for t in tiers_order:
        b = _as_float((tiers_raw.get(t) or {}).get("budget_s"))
        if b is not None and b > 0.0:
            tier_budgets[t] = b

    inst_meta: dict[str, dict[str, Any]] = {}
    insts_by_tier: dict[str, list[str]] = {t: [] for t in tiers_order}
    for inst in manifest.get("instances") or []:
        if not isinstance(inst, Mapping) or not inst.get("id"):
            continue
        iid = str(inst["id"])
        tier = str(inst.get("tier") or "T0")
        if tier not in insts_by_tier:
            tiers_order.append(tier)
            insts_by_tier[tier] = []
        insts_by_tier[tier].append(iid)
        br = inst.get("bracket_ref")
        bracket_norm = (
            [_as_int(br[0]), _as_int(br[1])]
            if isinstance(br, Sequence) and not isinstance(br, (str, bytes)) and len(br) == 2
            else None
        )
        inst_meta[iid] = {
            "id": iid,
            "tier": tier,
            "kind": inst.get("kind"),
            "d_ref": _as_int(inst.get("d_ref")),
            "bracket_ref": bracket_norm,
            "tags": dict(inst.get("tags") or {}),
            "in_manifest": True,
        }
    return tiers_order, tier_budgets, inst_meta, insts_by_tier


def _dedupe_and_filter(
    manifest: Mapping[str, Any],
    records: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    """
    Deduplicate records (last record per `run_id` wins) and separate manifest mismatches.

    Returns:
        `(valid_records, excluded_records, n_duplicates_overwritten)`.
    """
    by_id: dict[str, dict[str, Any]] = {}
    n_dups = 0
    for idx, rec in enumerate(records):
        if not isinstance(rec, Mapping):
            continue
        rid = rec.get("run_id")
        key = str(rid) if rid else f"__anon_{idx}"
        if key in by_id:
            n_dups += 1
        by_id[key] = dict(rec)

    m_id = manifest.get("id")
    m_ver = manifest.get("version")
    valid: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for rec in by_id.values():
        rm = rec.get("manifest") or {}
        rm_id = rm.get("id")
        rm_ver = rm.get("version")
        mismatch_reasons: list[str] = []
        if m_id is not None and rm_id is not None and str(rm_id) != str(m_id):
            mismatch_reasons.append(f"manifest.id={rm_id!r} != {m_id!r}")
        if m_ver is not None and rm_ver is not None and str(rm_ver) != str(m_ver):
            mismatch_reasons.append(f"manifest.version={rm_ver!r} != {m_ver!r}")
        if mismatch_reasons:
            excluded.append(
                {
                    "run_id": rec.get("run_id"),
                    "entry": _entry_name(rec),
                    "instance": (rec.get("instance") or {}).get("id"),
                    "reason": "; ".join(mismatch_reasons),
                }
            )
            continue
        valid.append(rec)
    return valid, excluded, n_dups


def _resolve_instance(
    rec: Mapping[str, Any],
    inst_meta: dict[str, dict[str, Any]],
    insts_by_tier: dict[str, list[str]],
    tiers_order: list[str],
) -> tuple[str, str]:
    """Return `(instance_id, tier)` using manifest metadata first, falling back to `record.instance`."""
    ri = rec.get("instance") or {}
    iid = str(ri.get("id") or "unknown")
    if iid in inst_meta:
        return iid, str(inst_meta[iid]["tier"])
    tier = str(ri.get("tier") or "T0")
    if tier not in insts_by_tier:
        tiers_order.append(tier)
        insts_by_tier[tier] = []
    insts_by_tier[tier].append(iid)
    br = ri.get("bracket_ref")
    bracket_norm = (
        [_as_int(br[0]), _as_int(br[1])]
        if isinstance(br, Sequence) and not isinstance(br, (str, bytes)) and len(br) == 2
        else None
    )
    inst_meta[iid] = {
        "id": iid,
        "tier": tier,
        "kind": ri.get("kind"),
        "d_ref": _as_int(ri.get("d_ref")),
        "bracket_ref": bracket_norm,
        "tags": {},
        "in_manifest": False,
    }
    return iid, tier


def _aggregate_cells(
    main_recs: Sequence[dict[str, Any]],
    inst_meta: dict[str, dict[str, Any]],
    insts_by_tier: dict[str, list[str]],
    tiers_order: list[str],
    tier_budgets: Mapping[str, float],
    baseline: str,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, list[str]], list[str]]:
    """Aggregate main-set records per `(entry, instance_id)`."""
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    entries_per_tier_set: dict[str, set[str]] = {t: set() for t in tiers_order}
    all_entries_set: set[str] = set()

    for rec in main_recs:
        iid, tier = _resolve_instance(rec, inst_meta, insts_by_tier, tiers_order)
        entry = _entry_name(rec)
        grouped.setdefault((entry, iid), []).append(rec)
        entries_per_tier_set.setdefault(tier, set()).add(entry)
        all_entries_set.add(entry)

    cells: dict[tuple[str, str], dict[str, Any]] = {}
    for (entry, iid), runs in grouped.items():
        tier = str(inst_meta[iid]["tier"])
        n_runs = len(runs)
        par2_vals = [_record_par2(r, tier_budgets, tier) for r in runs]
        par10_vals = [_record_par10(r, tier_budgets, tier) for r in runs]
        budgets = [_record_budget(r, tier_budgets, tier) for r in runs]
        solved_walls = [
            float((r.get("result") or {})["wall_sec"])
            for r in runs
            if _is_solved(r) and _as_float((r.get("result") or {}).get("wall_sec")) is not None
        ]
        n_solved = len(solved_walls)
        solved_by_med = (2 * n_solved) > n_runs
        med_par2 = float(np.median(par2_vals))
        med_par10 = float(np.median(par10_vals))
        min_par2 = float(min(par2_vals))
        max_par2 = float(max(par2_vals))
        t_min = float(min(solved_walls)) if solved_walls else None
        t_max = float(max(solved_walls)) if solved_walls else None
        t_med_solved = float(np.median(solved_walls)) if solved_walls else None

        lbs: list[int] = []
        ubs: list[int] = []
        status_counts: dict[str, int] = {}
        seeds: list[int] = []
        for r in runs:
            res = r.get("result") or {}
            st = str(res.get("status") or "unknown")
            status_counts[st] = status_counts.get(st, 0) + 1
            sd = _as_int((r.get("options") or {}).get("seed"))
            if sd is not None:
                seeds.append(sd)
            if st in VALID_BRACKET_STATUSES:
                lb = _as_int(res.get("lb"))
                ub = _as_int(res.get("ub"))
                if lb is not None:
                    lbs.append(lb)
                if ub is not None:
                    ubs.append(ub)

        best_lb = max(lbs) if lbs else None
        best_ub = min(ubs) if ubs else None
        best_gap = (best_ub - best_lb) if (best_lb is not None and best_ub is not None) else None

        cells[(entry, iid)] = {
            "entry": entry,
            "instance": iid,
            "tier": tier,
            "n_runs": n_runs,
            "n_solved": n_solved,
            "solved_by_median": solved_by_med,
            "budget_s": float(max(budgets)) if budgets else float(tier_budgets.get(tier, 0.0)),
            "med_par2": med_par2,
            "med_par10": med_par10,
            "min_par2": min_par2,
            "max_par2": max_par2,
            "t_med_solved": t_med_solved,
            "t_min": t_min,
            "t_max": t_max,
            "best_lb": best_lb,
            "best_ub": best_ub,
            "best_gap": best_gap,
            "status_counts": status_counts,
            "seeds": sorted(set(seeds)),
        }

    entries_by_tier = {t: _sort_entries(list(entries_per_tier_set.get(t, set())), baseline) for t in tiers_order}
    all_entries = _sort_entries(list(all_entries_set), baseline)
    return cells, entries_by_tier, all_entries


def _bootstrap_geomean_ci(
    ratios: Sequence[float],
    n_boot: int,
    seed: int,
) -> tuple[float | None, float | None, float | None, float | None, float | None]:
    """Compute `(geomean, ci_low, ci_high, min_r, max_r)` via bootstrap resampling over instances."""
    vals = [float(r) for r in ratios if r > 0.0 and math.isfinite(r)]
    n = len(vals)
    if n == 0:
        return None, None, None, None, None
    arr = np.asarray(vals, dtype=np.float64)
    log_arr = np.log(arr)
    gm = float(np.exp(np.mean(log_arr)))
    r_min = float(np.min(arr))
    r_max = float(np.max(arr))
    if n == 1 or n_boot <= 0:
        return gm, gm, gm, r_min, r_max
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot_gm = np.exp(np.mean(log_arr[idx], axis=1))
    ci_lo, ci_hi = np.percentile(boot_gm, [2.5, 97.5])
    return gm, float(ci_lo), float(ci_hi), r_min, r_max


def _classify_wtl(
    base_cell: Mapping[str, Any],
    entry_cell: Mapping[str, Any],
) -> str:
    """
    Classify an instance where both baseline and entry solved by median as `win`, `tie`, or `loss`.

    A tie holds when `max(t_b, t_e) <= 1.10 * min(t_b, t_e)` (within +-10%) OR when the seed
    min-max ranges `[min_par2, max_par2]` overlap (equal to solved `[t_min, t_max]` when all seeds
    are solved, and conservative when some seeds timed out).
    """
    tb = float(base_cell["med_par2"])
    te = float(entry_cell["med_par2"])
    lo_b, hi_b = float(base_cell["min_par2"]), float(base_cell["max_par2"])
    lo_e, hi_e = float(entry_cell["min_par2"]), float(entry_cell["max_par2"])

    t_lo = min(tb, te)
    t_hi = max(tb, te)
    within_10 = (t_hi <= 1.10 * t_lo) if t_lo > 0.0 else (t_hi == 0.0)
    overlap = (lo_e <= hi_b) and (lo_b <= hi_e)
    if within_10 or overlap:
        return "tie"
    return "win" if te < tb else "loss"


def _compute_vbs_for_instances(
    inst_ids: Sequence[str],
    entries: Sequence[str],
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    tier_budgets: Mapping[str, float],
    inst_meta: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, dict[str, Any]], dict[str, float | None], float | None]:
    """
    Compute per-instance VBS, per-entry marginal VBS contribution, and tier VBS mean PAR-2.

    For each instance `i` that has data from at least one entry:
    - `VBS(i) = min_e med_par2(e, i)`
    - `VBS_without(e, i) = min_{o != e} med_par2(o, i)` if another entry has data, else `2 * budget_s(i)`
    - `marginal_vbs(e) = mean_i VBS_without(e, i) - mean_i VBS(i)`
    """
    per_inst: dict[str, dict[str, Any]] = {}
    scored_insts: list[str] = []

    for iid in inst_ids:
        avail = [(e, cells[(e, iid)]) for e in entries if (e, iid) in cells]
        if not avail:
            per_inst[iid] = {
                "has_data": False,
                "solved_by_median": False,
                "med_par2": None,
                "best_entries": [],
                "best_lb": None,
                "best_ub": None,
                "best_gap": None,
                "marginal_by_entry": {},
            }
            continue
        scored_insts.append(iid)
        min_par2 = min(float(c["med_par2"]) for _, c in avail)
        solved_any = any(bool(c["solved_by_median"]) for _, c in avail)
        tol = max(1e-12, 1e-12 * abs(min_par2))
        best_entries = [e for e, c in avail if abs(float(c["med_par2"]) - min_par2) <= tol]

        lbs = [int(c["best_lb"]) for _, c in avail if c["best_lb"] is not None]
        ubs = [int(c["best_ub"]) for _, c in avail if c["best_ub"] is not None]
        vbs_lb = max(lbs) if lbs else None
        vbs_ub = min(ubs) if ubs else None
        vbs_gap = (vbs_ub - vbs_lb) if (vbs_lb is not None and vbs_ub is not None) else None

        tier = str(inst_meta[iid]["tier"])
        budget_i = max(float(c["budget_s"]) for _, c in avail)
        if budget_i <= 0.0:
            budget_i = float(tier_budgets.get(tier, 0.0))

        marg_i: dict[str, float] = {}
        for e in entries:
            others = [float(c["med_par2"]) for oe, c in avail if oe != e]
            without_e = min(others) if others else (2.0 * budget_i)
            marg_i[e] = float(without_e - min_par2)

        per_inst[iid] = {
            "has_data": True,
            "solved_by_median": solved_any,
            "med_par2": float(min_par2),
            "best_entries": best_entries,
            "best_lb": vbs_lb,
            "best_ub": vbs_ub,
            "best_gap": vbs_gap,
            "budget_s": budget_i,
            "marginal_by_entry": marg_i,
        }

    if not scored_insts:
        return per_inst, {e: None for e in entries}, None

    vbs_par2 = float(np.mean([float(per_inst[i]["med_par2"]) for i in scored_insts]))
    marg_by_entry: dict[str, float | None] = {}
    for e in entries:
        deltas = [float(per_inst[i]["marginal_by_entry"][e]) for i in scored_insts]
        marg_by_entry[e] = float(np.mean(deltas))
    return per_inst, marg_by_entry, vbs_par2


def _compute_table_a(
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    main_recs: Sequence[dict[str, Any]],
    tiers_order: Sequence[str],
    insts_by_tier: Mapping[str, list[str]],
    inst_meta: Mapping[str, Mapping[str, Any]],
    entries_by_tier: Mapping[str, list[str]],
    all_entries: Sequence[str],
    baseline: str,
    tier_budgets: Mapping[str, float],
    n_boot: int,
    boot_seed: int,
) -> dict[str, dict[str, Any]]:
    """Compute Table A statistics per tier (`T0`..`T3`) and combined (`all`), including a VBS row."""
    table_a: dict[str, dict[str, Any]] = {}
    scopes: list[tuple[str, list[str], list[str]]] = [
        (t, list(insts_by_tier.get(t, [])), list(entries_by_tier.get(t, [])))
        for t in tiers_order
    ]
    all_insts = [iid for t in tiers_order for iid in insts_by_tier.get(t, [])]
    scopes.append(("all", all_insts, list(all_entries)))

    for scope_name, scope_insts, scope_entries in scopes:
        scope_inst_set = set(scope_insts)
        n_manifest = sum(1 for i in scope_insts if inst_meta[i].get("in_manifest", True))
        n_total_insts = len(scope_insts)
        vbs_per_inst, marg_by_entry, vbs_par2 = _compute_vbs_for_instances(
            scope_insts, scope_entries, cells, tier_budgets, inst_meta
        )

        # Instances solved by ALL entries present in this scope.
        if scope_entries:
            all_solved_insts = [
                iid
                for iid in scope_insts
                if all(
                    (e, iid) in cells and bool(cells[(e, iid)]["solved_by_median"])
                    for e in scope_entries
                )
            ]
        else:
            all_solved_insts = []

        rows: list[dict[str, Any]] = []
        for entry in scope_entries:
            e_insts = [iid for iid in scope_insts if (entry, iid) in cells]
            runs_cert = sum(int(cells[(entry, i)]["n_solved"]) for i in e_insts)
            runs_tot = sum(int(cells[(entry, i)]["n_runs"]) for i in e_insts)
            inst_cert = sum(1 for i in e_insts if cells[(entry, i)]["solved_by_median"])
            inst_with_data = len(e_insts)
            inst_missing = n_total_insts - inst_with_data

            par2_mean = (
                float(np.mean([cells[(entry, i)]["med_par2"] for i in e_insts]))
                if e_insts
                else None
            )
            par10_mean = (
                float(np.mean([cells[(entry, i)]["med_par10"] for i in e_insts]))
                if e_insts
                else None
            )

            common_med = (
                float(np.median([cells[(entry, i)]["med_par2"] for i in all_solved_insts]))
                if all_solved_insts
                else None
            )

            # Pairwise vs baseline
            common_base = [i for i in e_insts if (baseline, i) in cells]
            par2_ratio = None
            if common_base:
                b_p2 = float(np.mean([cells[(baseline, i)]["med_par2"] for i in common_base]))
                e_p2 = float(np.mean([cells[(entry, i)]["med_par2"] for i in common_base]))
                if e_p2 > 0.0:
                    par2_ratio = b_p2 / e_p2

            both_solved = [
                i
                for i in common_base
                if cells[(baseline, i)]["solved_by_median"] and cells[(entry, i)]["solved_by_median"]
            ]
            only_entry = sum(
                1
                for i in common_base
                if cells[(entry, i)]["solved_by_median"] and not cells[(baseline, i)]["solved_by_median"]
            )
            only_base = sum(
                1
                for i in common_base
                if cells[(baseline, i)]["solved_by_median"] and not cells[(entry, i)]["solved_by_median"]
            )

            ratios = [
                float(cells[(baseline, i)]["med_par2"]) / float(cells[(entry, i)]["med_par2"])
                for i in both_solved
                if float(cells[(baseline, i)]["med_par2"]) > 0.0
                and float(cells[(entry, i)]["med_par2"]) > 0.0
            ]
            if entry == baseline and both_solved:
                gm, ci_lo, ci_hi, r_min, r_max = 1.0, 1.0, 1.0, 1.0, 1.0
            else:
                gm, ci_lo, ci_hi, r_min, r_max = _bootstrap_geomean_ci(
                    ratios, n_boot, _stable_seed(boot_seed, scope_name, entry)
                )

            wins = ties = losses = 0
            for i in both_solved:
                if entry == baseline:
                    ties += 1
                else:
                    outcome = _classify_wtl(cells[(baseline, i)], cells[(entry, i)])
                    if outcome == "win":
                        wins += 1
                    elif outcome == "tie":
                        ties += 1
                    else:
                        losses += 1

            rows.append(
                {
                    "entry": entry,
                    "is_baseline": entry == baseline,
                    "is_vbs": False,
                    "runs_certified": runs_cert,
                    "runs_total": runs_tot,
                    "instances_certified": inst_cert,
                    "instances_with_data": inst_with_data,
                    "instances_missing": inst_missing,
                    "instances_manifest": n_manifest,
                    "par2_mean": par2_mean,
                    "par10_mean": par10_mean,
                    "par2_n": inst_with_data,
                    "par2_ratio_vs_baseline": par2_ratio,
                    "par2_ratio_n": len(common_base),
                    "common_all_solved_median_s": common_med,
                    "common_all_solved_n": len(all_solved_insts),
                    "speedup_geomean": gm,
                    "speedup_ci_low": ci_lo,
                    "speedup_ci_high": ci_hi,
                    "speedup_min": r_min,
                    "speedup_max": r_max,
                    "speedup_n": len(ratios),
                    "wins": wins if common_base else None,
                    "ties": ties if common_base else None,
                    "losses": losses if common_base else None,
                    "only_entry": only_entry if common_base else None,
                    "only_baseline": only_base if common_base else None,
                    "marginal_vbs": marg_by_entry.get(entry),
                }
            )

        # VBS row for this scope
        vbs_row = None
        if scope_entries:
            # Count distinct (instance, seed) pairs certified by at least one entry in scope.
            pair_solved: dict[tuple[str, Any], bool] = {}
            for r in main_recs:
                iid = str((r.get("instance") or {}).get("id") or "")
                if iid not in scope_inst_set:
                    continue
                sd = (r.get("options") or {}).get("seed")
                pk = (iid, sd)
                pair_solved[pk] = pair_solved.get(pk, False) or _is_solved(r)
            vbs_runs_tot = len(pair_solved)
            vbs_runs_cert = sum(1 for ok in pair_solved.values() if ok)

            vbs_insts = [i for i in scope_insts if vbs_per_inst[i]["has_data"]]
            vbs_inst_cert = sum(1 for i in vbs_insts if vbs_per_inst[i]["solved_by_median"])
            vbs_common_med = (
                float(np.median([float(vbs_per_inst[i]["med_par2"]) for i in all_solved_insts]))
                if all_solved_insts
                else None
            )
            base_in_scope = [i for i in vbs_insts if (baseline, i) in cells]
            vbs_par2_ratio = None
            if base_in_scope:
                b_p2 = float(np.mean([cells[(baseline, i)]["med_par2"] for i in base_in_scope]))
                v_p2 = float(np.mean([float(vbs_per_inst[i]["med_par2"]) for i in base_in_scope]))
                if v_p2 > 0.0:
                    vbs_par2_ratio = b_p2 / v_p2
            vbs_both = [
                i
                for i in base_in_scope
                if cells[(baseline, i)]["solved_by_median"] and vbs_per_inst[i]["solved_by_median"]
            ]
            vbs_ratios = [
                float(cells[(baseline, i)]["med_par2"]) / float(vbs_per_inst[i]["med_par2"])
                for i in vbs_both
                if float(cells[(baseline, i)]["med_par2"]) > 0.0
                and float(vbs_per_inst[i]["med_par2"]) > 0.0
            ]
            v_gm, v_lo, v_hi, v_min, v_max = _bootstrap_geomean_ci(
                vbs_ratios, n_boot, _stable_seed(boot_seed, scope_name, "__VBS__")
            )
            v_wins = v_ties = v_losses = 0
            for i in vbs_both:
                best_e = vbs_per_inst[i]["best_entries"][0]
                if best_e == baseline:
                    v_ties += 1
                else:
                    oc = _classify_wtl(cells[(baseline, i)], cells[(best_e, i)])
                    if oc == "win":
                        v_wins += 1
                    elif oc == "tie":
                        v_ties += 1
                    else:
                        v_losses += 1
            v_only = sum(
                1
                for i in base_in_scope
                if vbs_per_inst[i]["solved_by_median"] and not cells[(baseline, i)]["solved_by_median"]
            )
            vbs_row = {
                "entry": "VBS",
                "is_baseline": False,
                "is_vbs": True,
                "runs_certified": vbs_runs_cert,
                "runs_total": vbs_runs_tot,
                "instances_certified": vbs_inst_cert,
                "instances_with_data": len(vbs_insts),
                "instances_missing": n_total_insts - len(vbs_insts),
                "instances_manifest": n_manifest,
                "par2_mean": vbs_par2,
                "par10_mean": None,
                "par2_n": len(vbs_insts),
                "par2_ratio_vs_baseline": vbs_par2_ratio,
                "par2_ratio_n": len(base_in_scope),
                "common_all_solved_median_s": vbs_common_med,
                "common_all_solved_n": len(all_solved_insts),
                "speedup_geomean": v_gm,
                "speedup_ci_low": v_lo,
                "speedup_ci_high": v_hi,
                "speedup_min": v_min,
                "speedup_max": v_max,
                "speedup_n": len(vbs_ratios),
                "wins": v_wins if base_in_scope else None,
                "ties": v_ties if base_in_scope else None,
                "losses": v_losses if base_in_scope else None,
                "only_entry": v_only if base_in_scope else None,
                "only_baseline": 0 if base_in_scope else None,
                "marginal_vbs": None,
            }

        obs_b = sorted(
            {
                float(cells[(e, i)]["budget_s"])
                for e in scope_entries
                for i in scope_insts
                if (e, i) in cells and cells[(e, i)].get("budget_s") is not None
            }
        )
        man_b = tier_budgets.get(scope_name)
        eff_b = max(obs_b) if (scope_name != "all" and obs_b) else man_b

        table_a[scope_name] = {
            "tier": scope_name,
            "budget_s": eff_b,
            "budget_manifest_s": man_b,
            "instances_manifest": n_manifest,
            "instances_total": n_total_insts,
            "entries": scope_entries,
            "all_solved_instances": all_solved_insts,
            "rows": rows,
            "vbs_row": vbs_row,
            "marginal_vbs": marg_by_entry,
            "has_data": bool(scope_entries),
        }

    return table_a


def _pick_best_entries_for_instance(
    entries: Sequence[str],
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    iid: str,
) -> list[str]:
    """Select the best entry (or tied entries) to bold in Table B for one instance."""
    avail = [(e, cells[(e, iid)]) for e in entries if (e, iid) in cells]
    if not avail:
        return []
    solved = [(e, c) for e, c in avail if c["solved_by_median"]]
    if solved:
        best_t = min(float(c["med_par2"]) for _, c in solved)
        tol = max(1e-12, 1e-12 * abs(best_t))
        return [e for e, c in solved if abs(float(c["med_par2"]) - best_t) <= tol]
    # Nobody solved by median: smallest valid bracket gap (tie-break: highest lb, lowest ub).
    with_gap = [(e, c) for e, c in avail if c["best_gap"] is not None]
    if with_gap:
        best_key = min((int(c["best_gap"]), -int(c["best_lb"]), int(c["best_ub"])) for _, c in with_gap)
        return [
            e
            for e, c in with_gap
            if (int(c["best_gap"]), -int(c["best_lb"]), int(c["best_ub"])) == best_key
        ]
    with_lb = [(e, c) for e, c in avail if c["best_lb"] is not None]
    if with_lb:
        max_lb = max(int(c["best_lb"]) for _, c in with_lb)
        return [e for e, c in with_lb if int(c["best_lb"]) == max_lb]
    return []


def _format_table_b_cell(cell: Mapping[str, Any] | None) -> str:
    """Format one `(entry, instance)` cell for Table B."""
    if cell is None:
        return "no data"
    n_runs = int(cell["n_runs"])
    n_solved = int(cell["n_solved"])
    if cell["solved_by_median"]:
        med_s = _fmt_time(cell["med_par2"])
        if n_runs == 1:
            return med_s
        lo_s = _fmt_time(cell["t_min"])
        hi_s = _fmt_time(cell["t_max"])
        base = f"{med_s} [{lo_s}, {hi_s}]"
        return base if n_solved == n_runs else f"{base} {n_solved}/{n_runs}"
    lb = cell["best_lb"]
    ub = cell["best_ub"]
    if lb is not None or ub is not None:
        lb_s = "-" if lb is None else str(int(lb))
        ub_s = "-" if ub is None else str(int(ub))
        return f"[{lb_s}, {ub_s}] {n_solved}/{n_runs}"
    counts = cell.get("status_counts") or {}
    if counts.get("wrong"):
        return f"wrong {n_solved}/{n_runs}"
    if counts.get("error"):
        return f"err {n_solved}/{n_runs}"
    return f"timeout {n_solved}/{n_runs}"


def _compute_table_b(
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    tiers_order: Sequence[str],
    insts_by_tier: Mapping[str, list[str]],
    inst_meta: Mapping[str, Mapping[str, Any]],
    entries_by_tier: Mapping[str, list[str]],
    tier_budgets: Mapping[str, float],
) -> dict[str, dict[str, Any]]:
    """Compute Table B per tier (manifest order): per-instance cells, best entry, VBS, and marginals."""
    table_b: dict[str, dict[str, Any]] = {}
    for tier in tiers_order:
        insts = list(insts_by_tier.get(tier, []))
        entries = list(entries_by_tier.get(tier, []))
        vbs_per_inst, marg_by_entry, vbs_par2 = _compute_vbs_for_instances(
            insts, entries, cells, tier_budgets, inst_meta
        )
        rows: list[dict[str, Any]] = []
        for iid in insts:
            meta = inst_meta[iid]
            best_entries = _pick_best_entries_for_instance(entries, cells, iid)
            entry_cells: dict[str, dict[str, Any] | None] = {}
            formatted: dict[str, str] = {}
            for e in entries:
                c = cells.get((e, iid))
                entry_cells[e] = dict(c) if c is not None else None
                formatted[e] = _format_table_b_cell(c)
            rows.append(
                {
                    "instance": iid,
                    "tier": tier,
                    "d_ref": meta.get("d_ref"),
                    "bracket_ref": meta.get("bracket_ref"),
                    "ref_display": _fmt_ref(meta.get("d_ref"), meta.get("bracket_ref")),
                    "cells": entry_cells,
                    "formatted": formatted,
                    "best_entries": best_entries,
                    "vbs": vbs_per_inst[iid],
                }
            )
        obs_b = sorted(
            {
                float(cells[(e, i)]["budget_s"])
                for e in entries
                for i in insts_by_tier.get(tier, [])
                if (e, i) in cells and cells[(e, i)].get("budget_s") is not None
            }
        )
        man_b = tier_budgets.get(tier)
        table_b[tier] = {
            "tier": tier,
            "budget_s": max(obs_b) if obs_b else man_b,
            "budget_manifest_s": man_b,
            "entries": entries,
            "rows": rows,
            "vbs_par2": vbs_par2,
            "marginal_vbs": marg_by_entry,
            "has_data": bool(entries),
        }
    return table_b


def _compute_table_c(
    valid_recs: Sequence[dict[str, Any]],
    excluded_recs: Sequence[dict[str, Any]],
    exclude_contaminated: bool,
) -> dict[str, Any]:
    """
    Compute Table C integrity audit across all deduplicated records matching the manifest.

    Checks:
    1. `wrong`: runs with `result.status == "wrong"` and their `verify.violations`.
    2. `lb_conflicts`: across all runs (all leagues and canaries), any run whose `result.lb` exceeds
       the minimum `verify.witness_weight` among runs of the same instance with
       `verify.witness_verified_numpy is True`.
    3. `unverified_certified`: runs with `result.status == "certified"` whose
       `verify.witness_verified_numpy` is not True.
    4. `contaminated`: runs with `env.contaminated` True.
    5. `errors`: runs with `result.status == "error"` (with `error.kind` and `error.msg`).
    6. `excluded`: records dropped before analysis (e.g. manifest id/version mismatch).
    """
    wrong_rows: list[dict[str, Any]] = []
    unverified_rows: list[dict[str, Any]] = []
    contaminated_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, Any]] = []

    # Collect minimum verified witness weight per instance across ALL valid records.
    min_witness: dict[str, dict[str, Any]] = {}
    for rec in valid_recs:
        iid = str((rec.get("instance") or {}).get("id") or "unknown")
        ver = rec.get("verify") or {}
        if ver.get("witness_verified_numpy") is True:
            w = _as_int(ver.get("witness_weight"))
            if w is not None:
                cur = min_witness.get(iid)
                if cur is None or w < int(cur["weight"]):
                    min_witness[iid] = {
                        "weight": w,
                        "run_id": rec.get("run_id"),
                        "entry": _entry_name(rec),
                    }

    lb_conflicts: list[dict[str, Any]] = []
    lb_conflicts_by_inst: dict[str, dict[str, Any]] = {}

    for rec in valid_recs:
        rid = rec.get("run_id")
        entry = _entry_name(rec)
        opt = rec.get("options") or {}
        league = str(opt.get("league") or "default")
        canary = opt.get("canary")
        seed = _as_int(opt.get("seed"))
        threads = _as_int(opt.get("threads"))
        iid = str((rec.get("instance") or {}).get("id") or "unknown")
        res = rec.get("result") or {}
        ver = rec.get("verify") or {}
        env_r = rec.get("env") or {}
        status = str(res.get("status") or "unknown")
        lb = _as_int(res.get("lb"))
        ub = _as_int(res.get("ub"))

        if status == "wrong":
            viols = ver.get("violations")
            v_list = [str(v) for v in viols] if isinstance(viols, Sequence) and not isinstance(viols, str) else []
            wrong_rows.append(
                {
                    "run_id": rid,
                    "entry": entry,
                    "league": league,
                    "canary": canary,
                    "instance": iid,
                    "seed": seed,
                    "lb": lb,
                    "ub": ub,
                    "violations": v_list,
                }
            )

        if iid in min_witness and lb is not None and lb > int(min_witness[iid]["weight"]):
            mw = min_witness[iid]
            lb_conflicts.append(
                {
                    "run_id": rid,
                    "entry": entry,
                    "league": league,
                    "canary": canary,
                    "instance": iid,
                    "seed": seed,
                    "status": status,
                    "lb": lb,
                    "min_verified_weight": int(mw["weight"]),
                    "witness_run_id": mw["run_id"],
                    "witness_entry": mw["entry"],
                }
            )
            inst_entry = lb_conflicts_by_inst.setdefault(
                iid,
                {
                    "instance": iid,
                    "max_lb": lb,
                    "min_verified_weight": int(mw["weight"]),
                    "witness_run_id": mw["run_id"],
                    "witness_entry": mw["entry"],
                    "n_conflicting_runs": 0,
                },
            )
            inst_entry["max_lb"] = max(int(inst_entry["max_lb"]), lb)
            inst_entry["n_conflicting_runs"] = int(inst_entry["n_conflicting_runs"]) + 1

        ub_ok = ver.get("ub_verified") if "ub_verified" in ver else ver.get("witness_verified_numpy")
        if status == "certified" and not res.get("k0") and ub_ok is not True:
            unverified_rows.append(
                {
                    "run_id": rid,
                    "entry": entry,
                    "league": league,
                    "canary": canary,
                    "instance": iid,
                    "seed": seed,
                    "ub": ub,
                    "ub_verified": ub_ok,
                    "witness_verified_numpy": ver.get("witness_verified_numpy"),
                    "unverified_reason": ver.get("unverified_reason"),
                }
            )

        if bool(env_r.get("contaminated")):
            contaminated_rows.append(
                {
                    "run_id": rid,
                    "entry": entry,
                    "league": league,
                    "canary": canary,
                    "instance": iid,
                    "seed": seed,
                    "threads": threads,
                    "wall_sec": _as_float(res.get("wall_sec")),
                    "ext_cpu_s": _as_float(env_r.get("ext_cpu_s")),
                    "steal_s": _as_float(env_r.get("steal_s")),
                    "excluded": bool(exclude_contaminated),
                }
            )

        if status == "error":
            err = res.get("error") or {}
            error_rows.append(
                {
                    "run_id": rid,
                    "entry": entry,
                    "league": league,
                    "canary": canary,
                    "instance": iid,
                    "seed": seed,
                    "kind": err.get("kind"),
                    "msg": err.get("msg"),
                }
            )

    return {
        "wrong": wrong_rows,
        "lb_conflicts": lb_conflicts,
        "lb_conflicts_by_instance": list(lb_conflicts_by_inst.values()),
        "unverified_certified": unverified_rows,
        "contaminated": contaminated_rows,
        "errors": error_rows,
        "excluded": list(excluded_recs),
        "exclude_contaminated": bool(exclude_contaminated),
        "clean": not (wrong_rows or lb_conflicts or unverified_rows or contaminated_rows),
    }


def _compute_table_d(
    valid_recs: Sequence[dict[str, Any]],
    envs_map: Mapping[str, dict[str, Any]],
) -> dict[str, Any]:
    """Compute Table D: environment provenance from env.json and engine versions from run records."""
    referenced_env_ids: set[str] = set()
    for rec in valid_recs:
        eid = (rec.get("env") or {}).get("env_id")
        if eid:
            referenced_env_ids.add(str(eid))

    env_rows: list[dict[str, Any]] = []
    for eid in sorted(envs_map):
        env = envs_map[eid]
        cpu = env.get("cpu") or {}
        mem = env.get("memory") or {}
        comp = env.get("compiler") or {}
        py = env.get("python") or {}
        qc = env.get("qec_cert") or {}
        dm = env.get("dist_m4ri") or {}
        harn = env.get("harness") or {}
        man = env.get("manifest") or {}
        inp = env.get("inputs") or {}
        sched = env.get("sched") or {}
        total_kb = _as_float(mem.get("total_kb"))
        mem_gib = (total_kb / (1024.0 * 1024.0)) if total_kb is not None else None
        env_rows.append(
            {
                "env_id": eid,
                "created": env.get("created"),
                "host": env.get("host"),
                "kernel": env.get("kernel"),
                "os": env.get("os"),
                "vm": env.get("vm"),
                "sched_autogroup": sched.get("autogroup"),
                "cpu": {
                    "model": cpu.get("model"),
                    "sockets": cpu.get("sockets"),
                    "cores_per_socket": cpu.get("cores_per_socket"),
                    "threads_per_core": cpu.get("threads_per_core"),
                    "cpus": cpu.get("cpus"),
                    "numa": cpu.get("numa"),
                    "caches": cpu.get("caches"),
                    "hypervisor": cpu.get("hypervisor"),
                    "flags_sha1": cpu.get("flags_sha1"),
                },
                "memory_gib": mem_gib,
                "memory_total_kb": total_kb,
                "compiler": {
                    "cc": comp.get("cc"),
                    "flags": comp.get("flags") if comp.get("flags") is not None else comp.get("makefile"),
                },
                "m4ri": env.get("m4ri"),
                "python": {
                    "version": py.get("version"),
                    "stim": py.get("stim"),
                    "numpy": py.get("numpy"),
                    "scipy": py.get("scipy"),
                    "matplotlib": py.get("matplotlib"),
                },
                "qec_cert": {
                    "version": qc.get("version"),
                    "commit": qc.get("commit"),
                    "dirty": qc.get("dirty"),
                    "lib_sha256": qc.get("lib_sha256"),
                },
                "dist_m4ri": {
                    "version": dm.get("version"),
                    "commit": dm.get("commit"),
                    "dirty": dm.get("dirty"),
                    "sha256": dm.get("sha256"),
                    "py_version": dm.get("py_version"),
                    "py_sha256": dm.get("py_sha256"),
                },
                "harness": {
                    "version": harn.get("version"),
                    "files": harn.get("files"),
                },
                "manifest": {
                    "id": man.get("id"),
                    "version": man.get("version"),
                    "sha256": man.get("sha256"),
                },
                "inputs": {
                    "checked": inp.get("checked"),
                    "files": inp.get("files"),
                    "bytes": inp.get("bytes"),
                    "sha256_ok": inp.get("sha256_ok"),
                    "page_cache_warmed": inp.get("page_cache_warmed"),
                },
            }
        )

    missing_env_ids = sorted(eid for eid in referenced_env_ids if eid not in envs_map)

    # Engine versions per entry across all valid records.
    eng_versions: dict[str, dict[tuple[Any, ...], dict[str, Any]]] = {}
    for rec in valid_recs:
        entry = _entry_name(rec)
        eng = rec.get("engine") or {}
        deps = eng.get("deps") or {}
        feats = eng.get("features") or {}
        vkey = (
            eng.get("name"),
            eng.get("interface"),
            eng.get("config"),
            eng.get("version"),
            eng.get("commit"),
            eng.get("dirty"),
            eng.get("sha256"),
        )
        bucket = eng_versions.setdefault(entry, {})
        if vkey not in bucket:
            bucket[vkey] = {
                "entry": entry,
                "name": eng.get("name"),
                "interface": eng.get("interface"),
                "config": eng.get("config"),
                "version": eng.get("version"),
                "commit": eng.get("commit"),
                "dirty": eng.get("dirty"),
                "sha256": eng.get("sha256"),
                "deps": dict(deps) if isinstance(deps, Mapping) else {},
                "features": dict(feats) if isinstance(feats, Mapping) else {},
                "n_runs": 0,
            }
        bucket[vkey]["n_runs"] = int(bucket[vkey]["n_runs"]) + 1

    def _entry_sort_key(ent: str) -> tuple[str, int]:
        if "@p" in ent:
            prefix, p_str = ent.rsplit("@p", 1)
            try:
                return (prefix, int(p_str))
            except ValueError:
                return (prefix, 10**9)
        return (ent, 10**9)

    engine_rows: list[dict[str, Any]] = []
    multi_version_entries: list[str] = []
    for entry in sorted(eng_versions, key=_entry_sort_key):
        variants = list(eng_versions[entry].values())
        conflict = len(variants) > 1
        if conflict:
            multi_version_entries.append(entry)
        for v in variants:
            v["version_conflict"] = conflict
            engine_rows.append(v)

    return {
        "envs": env_rows,
        "missing_env_ids": missing_env_ids,
        "engines": engine_rows,
        "multi_version_entries": multi_version_entries,
    }


def _compute_canary(canary_recs: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Compute canary drift statistics (wall_sec and t_stop over time)."""
    by_key: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for rec in canary_recs:
        entry = _entry_name(rec)
        iid = str((rec.get("instance") or {}).get("id") or "unknown")
        by_key.setdefault((entry, iid), []).append(rec)

    series_list: list[dict[str, Any]] = []
    for (entry, iid) in sorted(by_key):
        recs = sorted(by_key[(entry, iid)], key=lambda r: str((r.get("env") or {}).get("started_at") or ""))
        points: list[dict[str, Any]] = []
        walls: list[float] = []
        tstops: list[float] = []
        for idx, r in enumerate(recs):
            res = r.get("result") or {}
            ti = res.get("t_internal") or {}
            w = _as_float(res.get("wall_sec"))
            ts = _as_float(ti.get("t_stop"))
            if w is not None and _is_solved(r):
                walls.append(w)
            if ts is not None and _is_solved(r):
                tstops.append(ts)
            points.append(
                {
                    "index": idx + 1,
                    "canary": (r.get("options") or {}).get("canary"),
                    "started_at": (r.get("env") or {}).get("started_at"),
                    "status": res.get("status"),
                    "wall_sec": w,
                    "t_stop": ts,
                    "contaminated": bool((r.get("env") or {}).get("contaminated")),
                }
            )
        med_w = float(np.median(walls)) if walls else None
        min_w = float(min(walls)) if walls else None
        max_w = float(max(walls)) if walls else None
        spread = (100.0 * (max_w - min_w) / med_w) if (med_w and med_w > 0.0 and min_w is not None) else None
        med_ts = float(np.median(tstops)) if tstops else None
        series_list.append(
            {
                "entry": entry,
                "instance": iid,
                "n": len(recs),
                "n_certified": len(walls),
                "median_wall_s": med_w,
                "min_wall_s": min_w,
                "max_wall_s": max_w,
                "spread_pct": spread,
                "median_t_stop_s": med_ts,
                "points": points,
            }
        )
    return {"has_data": bool(series_list), "series": series_list}


def _compute_fig1_data(
    main_recs: Sequence[dict[str, Any]],
    manifest: Mapping[str, Any],
    tiers_order: Sequence[str],
    insts_by_tier: Mapping[str, list[str]],
    inst_meta: Mapping[str, Mapping[str, Any]],
    all_entries: Sequence[str],
    tier_budgets: Mapping[str, float],
) -> dict[str, Any]:
    """Compute Figure 1 (cactus/ECDF per tier + all tiers, with VBS curve and planned pairs)."""
    tiers_cfg = manifest.get("tiers") or {}
    panels: dict[str, dict[str, Any]] = {}
    scopes = list(tiers_order) + ["all"]

    for scope in scopes:
        if scope == "all":
            scope_insts = [i for t in tiers_order for i in insts_by_tier.get(t, [])]
            budget_s = max(tier_budgets.values()) if tier_budgets else None
            planned_pairs = 0
            for t in tiers_order:
                seeds_cfg = (tiers_cfg.get(t) or {}).get("seeds") or []
                planned_pairs += len(insts_by_tier.get(t, [])) * max(1, len(seeds_cfg))
        else:
            scope_insts = list(insts_by_tier.get(scope, []))
            budget_s = tier_budgets.get(scope)
            seeds_cfg = (tiers_cfg.get(scope) or {}).get("seeds") or []
            planned_pairs = len(scope_insts) * max(1, len(seeds_cfg)) if scope_insts else 0

        inst_set = set(scope_insts)
        scope_recs = [r for r in main_recs if str((r.get("instance") or {}).get("id") or "") in inst_set]
        obs_scope_b = [
            b
            for r in scope_recs
            if (b := _as_float((r.get("options") or {}).get("budget_s") or (r.get("result") or {}).get("cap_sec")))
            is not None
        ]
        if obs_scope_b:
            budget_s = max(obs_scope_b)

        curves: dict[str, list[float]] = {}
        pair_best: dict[tuple[str, Any], float] = {}
        seen_pairs: set[tuple[str, Any]] = set()
        for r in scope_recs:
            entry = _entry_name(r)
            iid = str((r.get("instance") or {}).get("id") or "")
            sd = (r.get("options") or {}).get("seed")
            pk = (iid, sd)
            seen_pairs.add(pk)
            if _is_solved(r):
                w = float((r.get("result") or {})["wall_sec"])
                curves.setdefault(entry, []).append(w)
                if pk not in pair_best or w < pair_best[pk]:
                    pair_best[pk] = w
            else:
                curves.setdefault(entry, [])

        for e in list(curves):
            curves[e].sort()
        vbs_times = sorted(pair_best.values())
        if planned_pairs == 0 and seen_pairs:
            planned_pairs = len(seen_pairs)

        panels[scope] = {
            "tier": scope,
            "budget_s": budget_s,
            "planned_pairs": planned_pairs,
            "observed_pairs": len(seen_pairs),
            "curves": {e: curves[e] for e in all_entries if e in curves},
            "vbs": vbs_times,
            "has_data": bool(scope_recs),
        }
    return {
        "tiers": list(tiers_order),
        "entries": list(all_entries),
        "panels": panels,
        "has_data": bool(main_recs),
    }


def _compute_fig2_data(
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    insts_by_tier: Mapping[str, list[str]],
    tiers_order: Sequence[str],
    all_entries: Sequence[str],
    baseline: str,
    tier_budgets: Mapping[str, float],
) -> dict[str, Any]:
    """Compute Figure 2 (baseline vs each non-baseline entry log-log scatter)."""
    others = [e for e in all_entries if e != baseline]
    all_insts = [i for t in tiers_order for i in insts_by_tier.get(t, [])]
    panels: list[dict[str, Any]] = []

    for entry in others:
        common = [i for i in all_insts if (baseline, i) in cells and (entry, i) in cells]
        if not common:
            panels.append(
                {"entry": entry, "baseline": baseline, "border": None, "points": [], "has_data": False}
            )
            continue
        max_budget = max(
            max(float(cells[(baseline, i)]["budget_s"]), float(cells[(entry, i)]["budget_s"]))
            for i in common
        )
        max_solved = 0.0
        for i in common:
            for e in (baseline, entry):
                if cells[(e, i)]["solved_by_median"]:
                    max_solved = max(max_solved, float(cells[(e, i)]["med_par2"]))
        border = max(2.0 * max_budget, 1.5 * max_solved, 1.0)
        points: list[dict[str, Any]] = []
        for i in common:
            cb = cells[(baseline, i)]
            ce = cells[(entry, i)]
            xb_solved = bool(cb["solved_by_median"])
            ye_solved = bool(ce["solved_by_median"])
            x_val = float(cb["med_par2"]) if xb_solved else border
            y_val = float(ce["med_par2"]) if ye_solved else border
            points.append(
                {
                    "instance": i,
                    "tier": cb["tier"],
                    "x_med_par2": float(cb["med_par2"]),
                    "y_med_par2": float(ce["med_par2"]),
                    "x_solved": xb_solved,
                    "y_solved": ye_solved,
                    "x_plot": x_val,
                    "y_plot": y_val,
                }
            )
        panels.append(
            {
                "entry": entry,
                "baseline": baseline,
                "border": border,
                "points": points,
                "has_data": True,
            }
        )
    return {
        "baseline": baseline,
        "tiers": list(tiers_order),
        "panels": panels,
        "has_data": any(p["has_data"] for p in panels),
    }


def _compute_fig3_data(
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    insts_by_tier: Mapping[str, list[str]],
    all_entries: Sequence[str],
    tier_budgets: Mapping[str, float],
) -> dict[str, Any]:
    """Compute Figure 3 (per-instance median PAR-2 bar + min-max whiskers on T1-T2)."""
    target_tiers = ("T1", "T2")
    inst_rows: list[dict[str, Any]] = []
    present_entries: set[str] = set()

    for tier in target_tiers:
        for iid in insts_by_tier.get(tier, []):
            bars: dict[str, dict[str, Any]] = {}
            for e in all_entries:
                c = cells.get((e, iid))
                if c is None:
                    continue
                present_entries.add(e)
                bars[e] = {
                    "med_par2": float(c["med_par2"]),
                    "min_par2": float(c["min_par2"]),
                    "max_par2": float(c["max_par2"]),
                    "solved_by_median": bool(c["solved_by_median"]),
                    "n_solved": int(c["n_solved"]),
                    "n_runs": int(c["n_runs"]),
                    "budget_s": float(c["budget_s"]),
                }
            if bars:
                inst_rows.append(
                    {
                        "instance": iid,
                        "tier": tier,
                        "budget_s": tier_budgets.get(tier),
                        "bars": bars,
                    }
                )
    ordered_entries = [e for e in all_entries if e in present_entries]
    return {
        "tiers": list(target_tiers),
        "entries": ordered_entries,
        "instances": inst_rows,
        "has_data": bool(inst_rows),
    }


def _compute_fig4_data(
    cells: Mapping[tuple[str, str], Mapping[str, Any]],
    insts_by_tier: Mapping[str, list[str]],
    inst_meta: Mapping[str, Mapping[str, Any]],
    all_entries: Sequence[str],
) -> dict[str, Any]:
    """Compute Figure 4 (T3 best brackets per entry with d_ref / bracket_ref markers)."""
    inst_rows: list[dict[str, Any]] = []
    present_entries: set[str] = set()

    for iid in insts_by_tier.get("T3", []):
        meta = inst_meta[iid]
        brackets: dict[str, dict[str, Any]] = {}
        for e in all_entries:
            c = cells.get((e, iid))
            if c is None:
                continue
            present_entries.add(e)
            brackets[e] = {
                "lb": c["best_lb"],
                "ub": c["best_ub"],
                "gap": c["best_gap"],
                "solved_by_median": bool(c["solved_by_median"]),
                "n_solved": int(c["n_solved"]),
                "n_runs": int(c["n_runs"]),
            }
        if brackets:
            inst_rows.append(
                {
                    "instance": iid,
                    "d_ref": meta.get("d_ref"),
                    "bracket_ref": meta.get("bracket_ref"),
                    "brackets": brackets,
                }
            )
    ordered_entries = [e for e in all_entries if e in present_entries]
    return {
        "entries": ordered_entries,
        "instances": inst_rows,
        "has_data": bool(inst_rows),
    }


def _is_scaling_complete(rec: Mapping[str, Any], config: str) -> bool:
    """Return True if a binary record is a complete scaling run (not killed, valid status & wall)."""
    res = rec.get("result") or {}
    if bool(res.get("killed")):
        return False
    w = _as_float(res.get("wall_sec"))
    if w is None or w <= 0.0:
        return False
    status = str(res.get("status") or "")
    if config == "cc_fixed":
        return status in ("certified", "bracket")
    return status == "certified"


def _compute_fig5_data(
    valid_recs: Sequence[dict[str, Any]],
    manifest: Mapping[str, Any],
    exclude_contaminated: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """
    Compute Figure 5 data: binary (`dist_m4ri`) thread scaling for `default` (method=3) and
    `cc_fixed` (method=2), combining league `scaling` records with league `default` records at
    p = 64 and p = 128 when `scaling` does not already cover that `(instance, default, p)`.
    """
    # Filter to binary, non-canary (and non-contaminated if excluded).
    scaling_recs: list[dict[str, Any]] = []
    default_p64_128: list[dict[str, Any]] = []

    for rec in valid_recs:
        eng = rec.get("engine") or {}
        opt = rec.get("options") or {}
        if eng.get("name") != "dist_m4ri":
            continue
        if opt.get("canary") is not None:
            continue
        if exclude_contaminated and bool((rec.get("env") or {}).get("contaminated")):
            continue
        league = str(opt.get("league") or "default")
        cfg = str(eng.get("config") or "default")
        thr = _as_int(opt.get("threads"))
        if thr is None or thr <= 0:
            continue
        if league == "scaling" and cfg in ("default", "cc_fixed"):
            scaling_recs.append(rec)
        elif league == "default" and cfg == "default" and thr in (64, 128):
            default_p64_128.append(rec)

    # Identify scaling instances: any instance present in league "scaling", plus manifest
    # scaling_instances that have >= 2 distinct thread counts across scaling + default p64/128.
    scaling_insts: list[str] = []
    seen_si: set[str] = set()
    manifest_si = [str(x) for x in (manifest.get("scaling_instances") or [])]

    for iid in manifest_si:
        if any((r.get("instance") or {}).get("id") == iid for r in scaling_recs):
            scaling_insts.append(iid)
            seen_si.add(iid)
    for r in scaling_recs:
        iid = str((r.get("instance") or {}).get("id") or "")
        if iid and iid not in seen_si:
            scaling_insts.append(iid)
            seen_si.add(iid)

    # Check manifest scaling_instances that only have default p64/128 runs.
    for iid in manifest_si:
        if iid in seen_si:
            continue
        ps = {
            _as_int((r.get("options") or {}).get("threads"))
            for r in default_p64_128
            if (r.get("instance") or {}).get("id") == iid and _is_scaling_complete(r, "default")
        }
        if len(ps) >= 2:
            scaling_insts.append(iid)
            seen_si.add(iid)

    # Fallback if no league=="scaling" and no manifest scaling_instances: use default instances with
    # both p=64 and p=128 complete.
    if not scaling_insts and default_p64_128:
        by_i_ps: dict[str, set[int]] = {}
        for r in default_p64_128:
            if not _is_scaling_complete(r, "default"):
                continue
            iid = str((r.get("instance") or {}).get("id") or "")
            thr = _as_int((r.get("options") or {}).get("threads"))
            if iid and thr is not None:
                by_i_ps.setdefault(iid, set()).add(thr)
        for iid in sorted(by_i_ps):
            if len(by_i_ps[iid]) >= 2:
                scaling_insts.append(iid)
                seen_si.add(iid)

    # Group scaling-league complete runs first; supplement missing (iid, "default", p in {64,128}).
    grouped_scaling: dict[tuple[str, str, int], list[dict[str, Any]]] = {}
    used_recs: list[dict[str, Any]] = []
    for r in scaling_recs:
        iid = str((r.get("instance") or {}).get("id") or "")
        cfg = str((r.get("engine") or {}).get("config") or "default")
        thr = _as_int((r.get("options") or {}).get("threads"))
        if iid in seen_si and thr is not None and _is_scaling_complete(r, cfg):
            grouped_scaling.setdefault((iid, cfg, thr), []).append(r)
            used_recs.append(r)

    grouped_default_supp: dict[tuple[str, str, int], list[dict[str, Any]]] = {}
    for r in default_p64_128:
        iid = str((r.get("instance") or {}).get("id") or "")
        thr = _as_int((r.get("options") or {}).get("threads"))
        if (
            iid in seen_si
            and thr is not None
            and (iid, "default", thr) not in grouped_scaling
            and _is_scaling_complete(r, "default")
        ):
            grouped_default_supp.setdefault((iid, "default", thr), []).append(r)
            used_recs.append(r)

    combined = {**grouped_scaling, **grouped_default_supp}

    series_list: list[dict[str, Any]] = []
    all_threads: set[int] = set()
    for iid in scaling_insts:
        for cfg in ("default", "cc_fixed"):
            ps = sorted(p for (si, sc, p) in combined if si == iid and sc == cfg)
            if not ps:
                continue
            p0 = ps[0]
            t0_walls = [float((r.get("result") or {})["wall_sec"]) for r in combined[(iid, cfg, p0)]]
            t0 = float(np.median(t0_walls))
            pts: list[dict[str, Any]] = []
            for p in ps:
                all_threads.add(p)
                ws = [float((r.get("result") or {})["wall_sec"]) for r in combined[(iid, cfg, p)]]
                med_t = float(np.median(ws))
                sp = (p0 * t0 / med_t) if med_t > 0.0 else None
                eff = (sp / p) if (sp is not None and p > 0) else None
                src = "scaling" if (iid, cfg, p) in grouped_scaling else "default"
                pts.append(
                    {
                        "threads": p,
                        "n": len(ws),
                        "med_wall_s": med_t,
                        "min_wall_s": float(min(ws)),
                        "max_wall_s": float(max(ws)),
                        "speedup": sp,
                        "efficiency": eff,
                        "source_league": src,
                    }
                )
            series_list.append(
                {
                    "instance": iid,
                    "config": cfg,
                    "method": 2 if cfg == "cc_fixed" else 3,
                    "p0": p0,
                    "t0": t0,
                    "points": pts,
                }
            )

    fig5 = {
        "instances": [i for i in scaling_insts if any(s["instance"] == i for s in series_list)],
        "threads": sorted(all_threads),
        "series": series_list,
        "has_data": bool(series_list),
    }
    return fig5, used_recs


def _engine_version_summary(recs: Sequence[dict[str, Any]], *, group_by_config: bool = False) -> str:
    """Summarize engines and versions observed in `recs` for figure captions."""
    if not recs:
        return "none (no data)"
    seen: dict[str, list[str]] = {}
    for r in recs:
        eng = r.get("engine") or {}
        if group_by_config:
            label = f"{eng.get('name') or 'unknown'}:{eng.get('config') or 'default'}"
        else:
            label = _entry_name(r)
        ver = str(eng.get("version") or "?")
        com = str(eng.get("commit") or "")
        dirty = "+dirty" if eng.get("dirty") else ""
        tag = f"{ver} ({com[:7]}{dirty})" if com else f"{ver}{dirty}"
        vals = seen.setdefault(label, [])
        if tag not in vals:
            vals.append(tag)
    parts = [f"{e} [{', '.join(vs)}]" for e, vs in sorted(seen.items())]
    return "; ".join(parts)


def _date_range(recs: Sequence[dict[str, Any]]) -> tuple[str | None, str | None, str]:
    """Compute `(min_started_at, max_ended_at, display_str)` over a record sequence."""
    starts: list[str] = []
    ends: list[str] = []
    for r in recs:
        ev = r.get("env") or {}
        s = ev.get("started_at")
        e = ev.get("ended_at") or s
        if s:
            starts.append(str(s))
        if e:
            ends.append(str(e))
    if not starts and not ends:
        return None, None, "unknown"
    lo = min(starts) if starts else min(ends)
    hi = max(ends) if ends else max(starts)
    lo_d = lo[:10]
    hi_d = hi[:10]
    disp = lo_d if lo_d == hi_d else f"{lo_d} to {hi_d}"
    return lo, hi, disp


def _seeds_summary(recs: Sequence[dict[str, Any]], inst_meta: Mapping[str, Mapping[str, Any]]) -> str:
    """Summarize observed seeds per tier (or overall if uniform) for captions."""
    if not recs:
        return "none"
    by_tier: dict[str, set[int]] = {}
    for r in recs:
        iid = str((r.get("instance") or {}).get("id") or "")
        tier = str((inst_meta.get(iid) or {}).get("tier") or (r.get("instance") or {}).get("tier") or "?")
        sd = _as_int((r.get("options") or {}).get("seed"))
        if sd is not None:
            by_tier.setdefault(tier, set()).add(sd)
    if not by_tier:
        return "unknown"

    def _fmt_set(ss: set[int]) -> str:
        xs = sorted(ss)
        if len(xs) >= 2 and xs == list(range(xs[0], xs[-1] + 1)):
            return f"{xs[0]}..{xs[-1]}"
        return ",".join(str(x) for x in xs)

    if len(by_tier) == 1:
        t, ss = next(iter(by_tier.items()))
        return f"{t}: {_fmt_set(ss)}"
    return ", ".join(f"{t}={_fmt_set(by_tier[t])}" for t in sorted(by_tier))


def _threads_summary(recs: Sequence[dict[str, Any]]) -> str:
    """Summarize observed thread counts for captions."""
    ps = sorted({p for r in recs if (p := _as_int((r.get("options") or {}).get("threads"))) is not None})
    if not ps:
        return "none"
    return ", ".join(str(p) for p in ps)


def _budgets_summary(
    recs: Sequence[dict[str, Any]],
    tier_budgets: Mapping[str, float],
    tiers: Sequence[str] | None = None,
) -> str:
    """Summarize observed budgets per tier (falling back to manifest budgets when no records exist)."""
    by_tier: dict[str, set[float]] = {}
    for r in recs:
        t = str((r.get("instance") or {}).get("tier") or "")
        b = _as_float((r.get("options") or {}).get("budget_s") or (r.get("result") or {}).get("cap_sec"))
        if t and b is not None:
            by_tier.setdefault(t, set()).add(b)

    parts: list[str] = []
    if tiers:
        active_tiers = [t for t in tiers if t in by_tier] if by_tier else list(tiers)
        for t in active_tiers:
            man_b = tier_budgets.get(t)
            obs = sorted(by_tier.get(t) or ())
            if obs:
                obs_s = "/".join(f"{_fmt_time(x)}s" for x in obs)
                if man_b is not None and (len(obs) != 1 or abs(obs[0] - man_b) > 1e-6):
                    parts.append(f"{t}={obs_s} (manifest {_fmt_time(man_b)}s)")
                else:
                    parts.append(f"{t}={obs_s}")
            elif man_b is not None:
                parts.append(f"{t}={_fmt_time(man_b)}s")
    if parts:
        return ", ".join(parts)
    bs = sorted({b for r in recs if (b := _as_float((r.get("options") or {}).get("budget_s"))) is not None})
    if not bs:
        return "unknown"
    return ", ".join(f"{_fmt_time(b)}s" for b in bs)


def _host_cpu_summary(envs_map: Mapping[str, dict[str, Any]]) -> str:
    """Summarize host and CPU model from envs for captions."""
    if not envs_map:
        return "host/CPU unknown (no env.json)"
    parts: list[str] = []
    for env in envs_map.values():
        host = env.get("host") or "unknown-host"
        cpu = env.get("cpu") or {}
        model = cpu.get("model") or "unknown CPU"
        cpus = cpu.get("cpus")
        s = f"{host} ({model}, {cpus} CPUs)" if cpus else f"{host} ({model})"
        if s not in parts:
            parts.append(s)
    return "; ".join(parts)


def _build_caption(
    title: str,
    recs: Sequence[dict[str, Any]],
    inst_meta: Mapping[str, Mapping[str, Any]],
    tier_budgets: Mapping[str, float],
    envs_map: Mapping[str, dict[str, Any]],
    how_to_read: str,
    *,
    tiers: Sequence[str] | None = None,
    exclude_contaminated: bool = False,
    n_contaminated: int = 0,
    group_by_config: bool = False,
) -> str:
    """Assemble a complete self-contained caption satisfying REVIEW.md 3.6."""
    eng_s = _engine_version_summary(recs, group_by_config=group_by_config)
    bud_s = _budgets_summary(recs, tier_budgets, tiers)
    thr_s = _threads_summary(recs)
    seed_s = _seeds_summary(recs, inst_meta)
    _, _, date_s = _date_range(recs)
    host_s = _host_cpu_summary(envs_map)
    contam_s = (
        f"Contaminated runs excluded ({n_contaminated} flagged)."
        if exclude_contaminated
        else (
            f"Includes {n_contaminated} contaminated run(s) (see Table C)."
            if n_contaminated > 0
            else "No contaminated runs detected."
        )
    )
    return (
        f"{title} Engines: {eng_s}. Budgets: {bud_s}. Threads: p={thr_s}. "
        f"Seeds: {seed_s}. {WALL_TIME_SENTENCE} Run dates: {date_s}. "
        f"Host: {host_s}. {contam_s} How to read: {how_to_read}"
    )


def analyze(
    manifest: Mapping[str, Any],
    records: Sequence[dict[str, Any]],
    envs: Mapping[str, dict[str, Any]] | Sequence[dict[str, Any]] | None = None,
    *,
    baseline: str = DEFAULT_BASELINE,
    exclude_contaminated: bool = False,
    n_boot: int = 2000,
    boot_seed: int = 1,
) -> dict[str, Any]:
    """
    Analyze certification benchmark records against a manifest (pure function; zero I/O).

    Deduplicates records (last record per `run_id` wins), separates manifest mismatches, filters
    main-set records (`league == "default"`, `canary is None`, optionally excluding contaminated
    runs), and computes Tables A-D, Figures 1-5 data, optional canary drift, and all captions.

    Args:
        manifest: Parsed `dm4ri-bench-manifest/1` dict.
        records: Sequence of `dm4ri-bench-run/1` dicts in session/file order.
        envs: Optional mapping `{env_id: env_dict}` or sequence of `env.json` dicts.
        baseline: Baseline engine entry (`name:config@p<threads>`), default `dist_m4ri:default@p128`.
        exclude_contaminated: If True, exclude runs with `env.contaminated == True` from main tables
            and figures (still reported in Table C).
        n_boot: Number of bootstrap resamples for geometric-mean speedup 95% CIs.
        boot_seed: Fixed RNG seed for deterministic bootstrap CIs.

    Returns:
        JSON-serializable summary dict with keys `schema`, `meta`, `entry_codes`, `table_a`,
        `table_b`, `table_c`, `table_d`, `canary`, `figures`, and `captions`.
    """
    envs_map = _normalize_envs(envs)
    tiers_order, tier_budgets, inst_meta, insts_by_tier = _build_manifest_index(manifest)
    valid_recs, excluded_recs, n_dups = _dedupe_and_filter(manifest, records)

    # Partition deduplicated manifest-matching records.
    main_recs: list[dict[str, Any]] = []
    canary_recs: list[dict[str, Any]] = []
    other_league_counts: dict[str, int] = {}
    n_contam_main = 0
    n_contam_all = 0

    for rec in valid_recs:
        opt = rec.get("options") or {}
        league = str(opt.get("league") or "default")
        is_canary = opt.get("canary") is not None
        is_contam = bool((rec.get("env") or {}).get("contaminated"))
        if is_contam:
            n_contam_all += 1
        if is_canary:
            if not (exclude_contaminated and is_contam):
                canary_recs.append(rec)
            continue
        if league == "default":
            if is_contam:
                n_contam_main += 1
            if exclude_contaminated and is_contam:
                continue
            main_recs.append(rec)
        else:
            other_league_counts[league] = other_league_counts.get(league, 0) + 1

    cells, entries_by_tier, all_entries = _aggregate_cells(
        main_recs, inst_meta, insts_by_tier, tiers_order, tier_budgets, baseline
    )

    table_a = _compute_table_a(
        cells,
        main_recs,
        tiers_order,
        insts_by_tier,
        inst_meta,
        entries_by_tier,
        all_entries,
        baseline,
        tier_budgets,
        n_boot,
        boot_seed,
    )
    table_b = _compute_table_b(cells, tiers_order, insts_by_tier, inst_meta, entries_by_tier, tier_budgets)
    table_c = _compute_table_c(valid_recs, excluded_recs, exclude_contaminated)
    table_d = _compute_table_d(valid_recs, envs_map)
    canary_data = _compute_canary(canary_recs)

    fig1 = _compute_fig1_data(main_recs, manifest, tiers_order, insts_by_tier, inst_meta, all_entries, tier_budgets)
    fig2 = _compute_fig2_data(cells, insts_by_tier, tiers_order, all_entries, baseline, tier_budgets)
    fig3 = _compute_fig3_data(cells, insts_by_tier, all_entries, tier_budgets)
    fig4 = _compute_fig4_data(cells, insts_by_tier, inst_meta, all_entries)
    fig5, scaling_used = _compute_fig5_data(valid_recs, manifest, exclude_contaminated)

    # Short entry codes (E0, E1, ...) for compact <=120-col Markdown tables.
    entry_codes = {e: f"E{idx}" for idx, e in enumerate(all_entries)}

    # Subset records for accurate per-figure captions.
    t1_t2_ids = set(insts_by_tier.get("T1", [])) | set(insts_by_tier.get("T2", []))
    t3_ids = set(insts_by_tier.get("T3", []))
    t1_t2_recs = [r for r in main_recs if str((r.get("instance") or {}).get("id") or "") in t1_t2_ids]
    t3_recs = [r for r in main_recs if str((r.get("instance") or {}).get("id") or "") in t3_ids]
    n_contam_t1_t2 = sum(
        1 for r in valid_recs
        if (r.get("options") or {}).get("league", "default") == "default"
        and (r.get("options") or {}).get("canary") is None
        and str((r.get("instance") or {}).get("id") or "") in t1_t2_ids
        and bool((r.get("env") or {}).get("contaminated"))
    )
    n_contam_t3 = sum(
        1 for r in valid_recs
        if (r.get("options") or {}).get("league", "default") == "default"
        and (r.get("options") or {}).get("canary") is None
        and str((r.get("instance") or {}).get("id") or "") in t3_ids
        and bool((r.get("env") or {}).get("contaminated"))
    )
    n_contam_scaling = sum(1 for r in scaling_used if bool((r.get("env") or {}).get("contaminated")))

    captions = {
        "fig1": _build_caption(
            "Figure 1: Cactus (ECDF) plot of certified (instance, seed) pairs vs. wall time by tier.",
            main_recs,
            inst_meta,
            tier_budgets,
            envs_map,
            (
                "Higher and further left is better. Each step curve counts certified (instance, seed) pairs "
                "finished within wall time t (log scale); dashed black shows the per-pair Virtual Best Solver "
                "(VBS); dotted horizontal line marks planned pairs; vertical dashed line marks the tier budget."
            ),
            tiers=tiers_order,
            exclude_contaminated=exclude_contaminated,
            n_contaminated=n_contam_main,
        ),
        "fig2": _build_caption(
            f"Figure 2: Pairwise log-log scatter of median wall time vs. baseline ({baseline}).",
            main_recs,
            inst_meta,
            tier_budgets,
            envs_map,
            (
                "Points below the y=x diagonal mean the entry is faster than the baseline; dotted lines mark "
                "10x faster (y=x/10) and 10x slower (y=10x). Instances unsolved by median sit on the dashed "
                "2*T_max border (open markers). Color encodes tier."
            ),
            tiers=tiers_order,
            exclude_contaminated=exclude_contaminated,
            n_contaminated=n_contam_main,
        ),
        "fig3": _build_caption(
            "Figure 3: Per-instance median PAR-2 wall time with min-max seed whiskers on T1-T2.",
            t1_t2_recs,
            inst_meta,
            tier_budgets,
            envs_map,
            (
                "Shorter bars are better (log time axis). Solid bars show median certified wall time across "
                "seeds with whiskers spanning [min, max] PAR-2; hatched bars mark instances not certified by "
                "the median (scored at 2*T)."
            ),
            tiers=("T1", "T2"),
            exclude_contaminated=exclude_contaminated,
            n_contaminated=n_contam_t1_t2,
        ),
        "fig4": _build_caption(
            "Figure 4: Best distance brackets [lb, ub] at the cap on T3 instances.",
            t3_recs,
            inst_meta,
            tier_budgets,
            envs_map,
            (
                "Narrower brackets are better; filled diamonds mark certified exact distances (lb == ub), "
                "bars span best [lb, ub] over seeds, right arrows denote lb with no upper bound found, "
                "black vertical ticks mark d_ref, and gray bands mark bracket_ref."
            ),
            tiers=("T3",),
            exclude_contaminated=exclude_contaminated,
            n_contaminated=n_contam_t3,
        ),
        "fig5": _build_caption(
            "Figure 5: Thread scaling of dist_m4ri binary: wall time T(p) and speedup S(p).",
            scaling_used,
            inst_meta,
            tier_budgets,
            envs_map,
            (
                "Left: median wall time T(p) with [min, max] whiskers (lower is better). Right: speedup "
                "S(p) = p0*T(p0)/T(p) normalized to each series' smallest thread count p0 against the dotted "
                "ideal line S(p)=p (higher is better). Vertical gray line at p=64 marks the SMT boundary; "
                "solid = default (method=3), dashed = cc_fixed (method=2)."
            ),
            exclude_contaminated=exclude_contaminated,
            n_contaminated=n_contam_scaling,
            group_by_config=True,
        ),
        "fig6": _build_caption(
            "Figure 6: Canary drift check (wall time and t_stop over session run order).",
            canary_recs,
            inst_meta,
            tier_budgets,
            envs_map,
            (
                "Flat points within the shaded +-10% band around the median indicate stable machine "
                "performance over time; open markers show engine search time t_stop when recorded."
            ),
            exclude_contaminated=exclude_contaminated,
            n_contaminated=sum(1 for r in canary_recs if bool((r.get("env") or {}).get("contaminated"))),
        ),
    }

    date_lo, date_hi, date_disp = _date_range(valid_recs)
    result = {
        "schema": SUMMARY_SCHEMA,
        "meta": {
            "manifest_id": manifest.get("id"),
            "manifest_version": manifest.get("version"),
            "baseline": baseline,
            "exclude_contaminated": bool(exclude_contaminated),
            "n_boot": int(n_boot),
            "boot_seed": int(boot_seed),
            "tiers": tiers_order,
            "tier_budgets": tier_budgets,
            "records_input": len(records),
            "records_duplicates_overwritten": n_dups,
            "records_excluded": len(excluded_recs),
            "records_valid": len(valid_recs),
            "records_main": len(main_recs),
            "records_canary": len(canary_recs),
            "records_other_leagues": other_league_counts,
            "contaminated_total": n_contam_all,
            "contaminated_main": n_contam_main,
            "date_start": date_lo,
            "date_end": date_hi,
            "date_range": date_disp,
            "wall_time_note": WALL_TIME_SENTENCE,
        },
        "entry_codes": entry_codes,
        "entries": all_entries,
        "table_a": table_a,
        "table_b": table_b,
        "table_c": table_c,
        "table_d": table_d,
        "canary": canary_data,
        "figures": {
            "fig1": fig1,
            "fig2": fig2,
            "fig3": fig3,
            "fig4": fig4,
            "fig5": fig5,
            "fig6": canary_data,
        },
        "captions": captions,
    }
    return _sanitize_json(result)


# ======================================================================================================================
# Figure rendering (matplotlib Agg, temp MPLCONFIGDIR)
# ======================================================================================================================

def _get_plt() -> Any:
    """Lazily import matplotlib.pyplot with Agg backend and a writable temp MPLCONFIGDIR."""
    global _MPL_TMPDIR, _PLT_MOD
    if _PLT_MOD is not None:
        return _PLT_MOD
    cur_cfg = os.environ.get("MPLCONFIGDIR")
    if not cur_cfg or not os.access(cur_cfg, os.W_OK):
        _MPL_TMPDIR = tempfile.mkdtemp(prefix="dm4ri-mpl-")
        os.environ["MPLCONFIGDIR"] = _MPL_TMPDIR
        atexit.register(shutil.rmtree, _MPL_TMPDIR, True)
    import matplotlib  # noqa: PLC0415

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    _PLT_MOD = plt
    return plt


def _draw_no_data(ax: Any, title: str, msg: str = "no data") -> None:
    """Render a clean placeholder axis when no records exist for a panel."""
    ax.set_title(title, fontsize=9)
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11, color="0.45", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("0.8")


def _entry_color_map(entries: Sequence[str]) -> dict[str, str]:
    """Assign stable Tab10 colors to entries."""
    palette = [f"C{i % 10}" for i in range(max(1, len(entries)))]
    return {e: palette[i] for i, e in enumerate(entries)}


def _render_fig1(
    fig_data: Mapping[str, Any],
    out_path: Path,
    color_map: Mapping[str, str] | None = None,
) -> None:
    """Render Figure 1: 2x3 cactus/ECDF grid (T0-T3, all tiers, shared legend)."""
    plt = _get_plt()
    tiers = list(fig_data.get("tiers") or TIER_ORDER_DEFAULT)
    entries = list(fig_data.get("entries") or [])
    panels = fig_data.get("panels") or {}
    colors = dict(color_map) if color_map is not None else _entry_color_map(entries)
    for idx_e, e in enumerate(entries):
        colors.setdefault(e, f"C{idx_e % 10}")
    scopes = tiers[:4] + ["all"]

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6))
    axes_flat = list(axes.flat)
    legend_handles: dict[str, Any] = {}

    for idx, scope in enumerate(scopes):
        ax = axes_flat[idx]
        pdata = panels.get(scope) or {}
        budget = pdata.get("budget_s")
        title = (
            f"{scope} (budget {_fmt_time(budget)} s)"
            if scope != "all" and budget is not None
            else ("All tiers" if scope == "all" else scope)
        )
        if not pdata.get("has_data"):
            _draw_no_data(ax, title)
            continue

        all_ts: list[float] = []
        for ts in (pdata.get("curves") or {}).values():
            all_ts.extend(float(x) for x in ts)
        all_ts.extend(float(x) for x in (pdata.get("vbs") or []))
        min_t = max(min(all_ts) * 0.6, 1e-4) if all_ts else 0.01
        max_t = max(max(all_ts) * 1.2, float(budget or 1.0)) if all_ts else float(budget or 60.0)
        if max_t <= min_t:
            max_t = min_t * 10.0

        for e in entries:
            ts = (pdata.get("curves") or {}).get(e)
            if ts is None:
                continue
            xs = [min_t] + [max(float(t), 1e-4) for t in ts] + [max_t]
            ys = [0] + list(range(1, len(ts) + 1)) + [len(ts)]
            (line,) = ax.step(xs, ys, where="post", color=colors[e], lw=1.8, label=e)
            legend_handles.setdefault(e, line)

        vbs = pdata.get("vbs") or []
        if vbs:
            xs = [min_t] + [max(float(t), 1e-4) for t in vbs] + [max_t]
            ys = [0] + list(range(1, len(vbs) + 1)) + [len(vbs)]
            (vline,) = ax.step(xs, ys, where="post", color="k", ls="--", lw=2.0, label="VBS")
            legend_handles.setdefault("VBS", vline)

        planned = int(pdata.get("planned_pairs") or 0)
        if planned > 0:
            hline = ax.axhline(planned, color="0.55", ls=":", lw=1.1, label="planned pairs")
            legend_handles.setdefault("planned pairs", hline)
        if scope != "all" and budget is not None and budget > 0:
            bline = ax.axvline(float(budget), color="0.65", ls="-.", lw=1.0, label="tier budget T")
            legend_handles.setdefault("tier budget T", bline)

        ax.set_xscale("log")
        ax.set_xlim(min_t, max_t * 1.05)
        ymax = max(planned, len(vbs), max((len(v) for v in (pdata.get("curves") or {}).values()), default=0), 1)
        ax.set_ylim(0, ymax * 1.12)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("wall time t (s)", fontsize=8.5)
        ax.set_ylabel("certified (instance, seed) pairs", fontsize=8.5)
        ax.grid(True, which="both", alpha=0.25)

    # 6th axis: shared legend in canonical order (entries, then VBS and reference lines).
    leg_ax = axes_flat[5]
    leg_ax.axis("off")
    if legend_handles:
        ordered_keys = [e for e in entries if e in legend_handles]
        for rk in ("VBS", "planned pairs", "tier budget T"):
            if rk in legend_handles:
                ordered_keys.append(rk)
        leg_ax.legend(
            [legend_handles[k] for k in ordered_keys],
            ordered_keys,
            loc="center",
            frameon=True,
            fontsize=8.5,
            title="Entries & reference",
            title_fontsize=9.5,
        )
    else:
        _draw_no_data(leg_ax, "Legend", "no data")

    fig.suptitle("Figure 1: Certification Cactus (ECDF) by Tier", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _render_fig2(fig_data: Mapping[str, Any], out_path: Path) -> None:
    """Render Figure 2: log-log scatter of baseline vs each non-baseline entry."""
    plt = _get_plt()
    panels = list(fig_data.get("panels") or [])
    baseline = str(fig_data.get("baseline") or DEFAULT_BASELINE)
    tier_colors = {"T0": "C0", "T1": "C1", "T2": "C2", "T3": "C3"}

    if not panels or not fig_data.get("has_data"):
        fig, ax = plt.subplots(figsize=(6.0, 5.2))
        _draw_no_data(ax, f"Baseline ({baseline}) vs. other entries", "no data (baseline or peer missing)")
        fig.tight_layout()
        fig.savefig(out_path, dpi=130)
        plt.close(fig)
        return

    n = len(panels)
    ncols = min(n, 3)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.7 * nrows), squeeze=False)
    axes_flat = list(axes.flat)

    for idx, pdata in enumerate(panels):
        ax = axes_flat[idx]
        entry = str(pdata["entry"])
        if not pdata.get("has_data"):
            _draw_no_data(ax, f"{entry} vs {baseline}")
            continue
        border = float(pdata["border"])
        pts = pdata["points"]
        pos_vals = [
            max(float(v), 1e-4)
            for p in pts
            for v in (p["x_plot"], p["y_plot"])
            if v is not None
        ]
        lo = max(min(pos_vals) * 0.5, 1e-4) if pos_vals else 1e-2
        hi = border * 1.45

        # Reference lines
        diag = np.geomspace(lo, hi, 50)
        ax.plot(diag, diag, color="0.25", lw=1.1, label="y = x")
        ax.plot(diag, diag * 10.0, color="0.55", ls=":", lw=1.0, label="10x / 0.1x")
        ax.plot(diag, diag / 10.0, color="0.55", ls=":", lw=1.0)
        ax.axvline(border, color="0.65", ls="--", lw=0.9)
        ax.axhline(border, color="0.65", ls="--", lw=0.9)

        seen_tiers: set[str] = set()
        for p in pts:
            tier = str(p["tier"])
            col = tier_colors.get(tier, "C4")
            xp = max(float(p["x_plot"]), 1e-4)
            yp = max(float(p["y_plot"]), 1e-4)
            both_ok = bool(p["x_solved"] and p["y_solved"])
            lab = tier if tier not in seen_tiers else None
            seen_tiers.add(tier)
            if both_ok:
                ax.scatter([xp], [yp], color=col, s=36, alpha=0.88, edgecolors="k", linewidths=0.4, label=lab)
            else:
                ax.scatter(
                    [xp], [yp], facecolors="none", edgecolors=col, marker="^", s=52, linewidths=1.3, label=lab
                )

        ax.text(border * 0.92, lo * 1.25, "unsolved\n(2T_max)", ha="right", va="bottom", fontsize=7.5, color="0.45")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{baseline} median time (s)", fontsize=8.5)
        ax.set_ylabel(f"{entry} median time (s)", fontsize=8.5)
        ax.set_title(entry, fontsize=9.5)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7.5, loc="upper left")

    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.suptitle(f"Figure 2: Pairwise Wall Time vs. Baseline ({baseline})", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _render_fig3(
    fig_data: Mapping[str, Any],
    out_path: Path,
    color_map: Mapping[str, str] | None = None,
) -> None:
    """Render Figure 3: grouped horizontal median PAR-2 bars + min-max whiskers on T1-T2."""
    plt = _get_plt()
    from matplotlib.patches import Patch  # noqa: PLC0415

    insts = list(fig_data.get("instances") or [])
    entries = list(fig_data.get("entries") or [])

    if not insts or not entries:
        fig, ax = plt.subplots(figsize=(8.0, 4.0))
        _draw_no_data(ax, "Figure 3: T1-T2 Per-Instance Wall Time (PAR-2)")
        fig.tight_layout()
        fig.savefig(out_path, dpi=130)
        plt.close(fig)
        return

    n_inst = len(insts)
    n_ent = len(entries)
    colors = dict(color_map) if color_map is not None else _entry_color_map(entries)
    for idx_e, e in enumerate(entries):
        colors.setdefault(e, f"C{idx_e % 10}")
    bar_h = min(0.78 / max(n_ent, 1), 0.34)
    fig_h = max(4.6, 0.36 * n_inst * max(1.0, 0.55 * n_ent) + 2.2)
    fig, ax = plt.subplots(figsize=(11.0, min(fig_h, 24.0)))

    all_vals: list[float] = []
    y_centers = np.arange(n_inst, dtype=np.float64)
    offsets = (np.arange(n_ent, dtype=np.float64) - 0.5 * (n_ent - 1)) * bar_h

    hatched_seen = False

    for i, item in enumerate(insts):
        bars = item["bars"]
        for j, e in enumerate(entries):
            b = bars.get(e)
            if b is None:
                continue
            y = y_centers[i] + offsets[j]
            med = max(float(b["med_par2"]), 1e-4)
            lo = max(float(b["min_par2"]), 1e-4)
            hi = max(float(b["max_par2"]), 1e-4)
            all_vals.extend([lo, med, hi])
            solved = bool(b["solved_by_median"])
            if solved:
                ax.barh(y, med, height=bar_h * 0.9, color=colors[e], alpha=0.85)
            else:
                hatched_seen = True
                ax.barh(
                    y,
                    med,
                    height=bar_h * 0.9,
                    color=colors[e],
                    alpha=0.35,
                    edgecolor=colors[e],
                    hatch="//",
                )
            if int(b["n_runs"]) > 1 and hi > lo:
                err_lo = max(0.0, med - lo)
                err_hi = max(0.0, hi - med)
                ax.errorbar(
                    med, y, xerr=[[err_lo], [err_hi]], fmt="none", ecolor="k", capsize=2.5, lw=0.9
                )

    # Separator line between T1 and T2 when both present.
    t1_count = sum(1 for item in insts if item["tier"] == "T1")
    if 0 < t1_count < n_inst:
        ax.axhline(t1_count - 0.5, color="0.4", ls="--", lw=1.0)

    labels = [f"[{item['tier']}] {_short_inst_label(item['instance'], 44)}" for item in insts]
    ax.set_yticks(y_centers)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(n_inst - 0.5, -0.5)
    ax.set_xscale("log")
    if all_vals:
        ax.set_xlim(max(min(all_vals) * 0.5, 1e-4), max(all_vals) * 1.6)
    ax.set_xlabel("median wall time / PAR-2 (s) with [min, max] whiskers over seeds", fontsize=9)
    ax.grid(True, axis="x", which="both", alpha=0.25)

    handles = [Patch(facecolor=colors[e], alpha=0.85, label=e) for e in entries]
    if hatched_seen:
        handles.append(
            Patch(facecolor="0.75", alpha=0.5, edgecolor="0.3", hatch="//", label="unsolved by median (2T)")
        )
    ax.legend(
        handles=handles,
        fontsize=8,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=min(len(handles), 3),
        frameon=True,
    )
    fig.suptitle(
        "Figure 3: Per-Instance Certification Time on T1-T2 (hatched = unsolved by median, 2T)",
        fontsize=10.5,
        y=0.99,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _render_fig4(
    fig_data: Mapping[str, Any],
    out_path: Path,
    color_map: Mapping[str, str] | None = None,
) -> None:
    """Render Figure 4: T3 distance brackets [lb, ub] per entry with d_ref / bracket_ref markers."""
    plt = _get_plt()
    from matplotlib.lines import Line2D  # noqa: PLC0415
    from matplotlib.patches import Patch  # noqa: PLC0415

    insts = list(fig_data.get("instances") or [])
    entries = list(fig_data.get("entries") or [])

    if not insts or not entries:
        fig, ax = plt.subplots(figsize=(8.5, 4.0))
        _draw_no_data(ax, "Figure 4: T3 Best Distance Brackets at Cap")
        fig.tight_layout()
        fig.savefig(out_path, dpi=130)
        plt.close(fig)
        return

    n_inst = len(insts)
    n_ent = len(entries)
    colors = dict(color_map) if color_map is not None else _entry_color_map(entries)
    for idx_e, e in enumerate(entries):
        colors.setdefault(e, f"C{idx_e % 10}")
    step = min(0.72 / max(n_ent, 1), 0.25)
    offsets = (np.arange(n_ent, dtype=np.float64) - 0.5 * (n_ent - 1)) * step
    fig_h = max(4.2, 0.56 * n_inst * max(1.0, 0.45 * n_ent) + 2.0)
    fig, ax = plt.subplots(figsize=(10.5, min(fig_h, 20.0)))

    # Collect x extent to size open rays nicely.
    x_vals: list[float] = []
    for item in insts:
        if item.get("d_ref") is not None:
            x_vals.append(float(item["d_ref"]))
        br = item.get("bracket_ref")
        if br and len(br) == 2:
            for v in br:
                if v is not None:
                    x_vals.append(float(v))
        for b in item["brackets"].values():
            if b.get("lb") is not None:
                x_vals.append(float(b["lb"]))
            if b.get("ub") is not None:
                x_vals.append(float(b["ub"]))
    x_max = max(x_vals) if x_vals else 20.0
    x_min = min(x_vals) if x_vals else 0.0
    ray_len = max(2.0, 0.08 * max(x_max - x_min, 10.0))

    dref_seen = False
    bref_seen = False
    cert_seen = False

    for i, item in enumerate(insts):
        if i > 0:
            ax.axhline(i - 0.5, color="0.88", lw=0.8)
        d_ref = item.get("d_ref")
        b_ref = item.get("bracket_ref")
        if d_ref is not None:
            dref_seen = True
            ax.plot([float(d_ref), float(d_ref)], [i - 0.42, i + 0.42], color="k", lw=2.0, alpha=0.8)
        elif b_ref and len(b_ref) == 2 and b_ref[0] is not None and b_ref[1] is not None:
            bref_seen = True
            ax.fill_betweenx([i - 0.42, i + 0.42], float(b_ref[0]), float(b_ref[1]), color="0.75", alpha=0.35)

        for j, e in enumerate(entries):
            b = item["brackets"].get(e)
            if b is None:
                continue
            y = i + offsets[j]
            lb = b.get("lb")
            ub = b.get("ub")
            col = colors[e]
            if lb is not None and ub is not None:
                if lb == ub:
                    ax.scatter([float(lb)], [y], color=col, marker="D", s=48, zorder=4)
                    cert_seen = True
                else:
                    ax.plot([float(lb), float(ub)], [y, y], color=col, lw=3.0, solid_capstyle="butt")
                    ax.scatter([float(lb), float(ub)], [y, y], color=col, marker="|", s=70, lw=1.8, zorder=4)
            elif lb is not None:
                ax.annotate(
                    "",
                    xy=(float(lb) + ray_len, y),
                    xytext=(float(lb), y),
                    arrowprops={"arrowstyle": "->", "color": col, "lw": 2.0},
                )
                ax.plot([float(lb)], [y], color=col, marker="|", ms=8, mew=1.8)
                ax.text(float(lb) + ray_len * 1.08, y, f">={lb}", va="center", fontsize=7.5, color=col)
            elif ub is not None:
                x0 = max(0.0, float(ub) - ray_len)
                ax.annotate(
                    "",
                    xy=(x0, y),
                    xytext=(float(ub), y),
                    arrowprops={"arrowstyle": "->", "color": col, "lw": 2.0},
                )
                ax.plot([float(ub)], [y], color=col, marker="|", ms=8, mew=1.8)

    handles: list[Any] = [Line2D([0], [0], color=colors[e], lw=3.0, label=e) for e in entries]
    if cert_seen:
        handles.append(Line2D([0], [0], color="0.3", marker="D", ls="none", ms=6, label="certified (lb = ub)"))
    if dref_seen:
        handles.append(Line2D([0], [0], color="k", lw=2.0, label="d_ref"))
    if bref_seen:
        handles.append(Patch(facecolor="0.75", alpha=0.45, label="bracket_ref"))

    ax.set_yticks(range(n_inst))
    ax.set_yticklabels([_short_inst_label(item["instance"], 44) for item in insts], fontsize=8)
    ax.set_ylim(n_inst - 0.5, -0.5)
    ax.set_xlim(max(0.0, x_min - 2.0), x_max + ray_len * 2.2 + 2.0)
    ax.set_xlabel("distance bound [lb, ub]", fontsize=9)
    ax.grid(True, axis="x", alpha=0.25)
    ax.legend(
        handles=handles,
        fontsize=8,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=min(len(handles), 3),
        frameon=True,
    )
    fig.suptitle("Figure 4: T3 Best Distance Brackets per Entry (with d_ref / bracket_ref)", fontsize=10.5, y=0.99)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _render_fig5(fig_data: Mapping[str, Any], out_path: Path) -> None:
    """Render Figure 5: dist_m4ri binary thread scaling — T(p) and speedup S(p)."""
    plt = _get_plt()
    series_list = list(fig_data.get("series") or [])

    if not series_list:
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
        _draw_no_data(axes[0], "Wall time T(p)")
        _draw_no_data(axes[1], "Speedup S(p)")
        fig.suptitle("Figure 5: Thread Scaling of dist_m4ri", fontsize=11)
        fig.tight_layout()
        fig.savefig(out_path, dpi=130)
        plt.close(fig)
        return

    insts = list(fig_data.get("instances") or [])
    inst_colors = {iid: f"C{i % 10}" for i, iid in enumerate(insts)}
    all_ps = sorted({int(pt["threads"]) for s in series_list for pt in s["points"]})

    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(12.2, 5.0))
    for ser in series_list:
        iid = str(ser["instance"])
        cfg = str(ser["config"])
        p0 = int(ser["p0"])
        pts = ser["points"]
        ps = [int(p["threads"]) for p in pts]
        ms = [max(float(p["med_wall_s"]), 1e-4) for p in pts]
        lo = [max(0.0, float(p["med_wall_s"]) - float(p["min_wall_s"])) for p in pts]
        hi = [max(0.0, float(p["max_wall_s"]) - float(p["med_wall_s"])) for p in pts]
        sps = [float(p["speedup"]) for p in pts if p["speedup"] is not None]
        ps_sp = [int(p["threads"]) for p in pts if p["speedup"] is not None]

        col = inst_colors.get(iid, "C0")
        ls = "-" if cfg == "default" else "--"
        mk = "o" if cfg == "default" else "s"
        cfg_lbl = "m3 (default)" if cfg == "default" else "m2 (cc_fixed)"
        p0_note = "" if p0 == 1 else f", p0={p0}"
        label = f"{_short_inst_label(iid, 26)} [{cfg_lbl}{p0_note}]"

        ax_t.errorbar(ps, ms, yerr=[lo, hi], fmt=f"{mk}{ls}", color=col, capsize=3, lw=1.6, ms=5, label=label)
        if ps_sp:
            ax_s.plot(ps_sp, sps, f"{mk}{ls}", color=col, lw=1.6, ms=5, label=label)

    p_lo = min(all_ps)
    p_hi = max(all_ps)
    p_line_lo = p_lo if p_lo < p_hi else max(1, p_lo // 2)
    p_line_hi = p_hi if p_lo < p_hi else p_hi * 2
    ax_s.plot([p_line_lo, p_line_hi], [p_line_lo, p_line_hi], "k:", lw=1.2, label="ideal S(p)=p")

    ticks = [p for p in (1, 2, 4, 8, 16, 32, 64, 128, 256) if p_line_lo <= p <= p_line_hi]
    if not ticks:
        ticks = all_ps

    for ax in (ax_t, ax_s):
        ax.set_xscale("log", base=2)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(p) for p in ticks])
        if p_line_lo < 64 <= p_line_hi:
            ax.axvline(64, color="0.6", ls="-.", lw=1.0)
        ax.grid(True, which="both", alpha=0.25)
        ax.set_xlabel("threads p (p > 64 uses SMT siblings)", fontsize=8.5)

    ax_t.set_yscale("log")
    ax_t.set_ylabel("median wall time T(p) (s)", fontsize=8.5)
    ax_t.set_title("Wall time T(p) with [min, max] whiskers", fontsize=9.5)
    ax_t.legend(fontsize=7, loc="best")

    ax_s.set_yscale("log", base=2)
    ax_s.set_yticks(ticks)
    ax_s.set_yticklabels([str(p) for p in ticks])
    ax_s.set_ylabel("speedup S(p) = p0 * T(p0) / T(p)", fontsize=8.5)
    ax_s.set_title("Parallel speedup S(p) (dotted = ideal)", fontsize=9.5)
    ax_s.legend(fontsize=7, loc="upper left")

    fig.suptitle("Figure 5: Thread Scaling of dist_m4ri (method=3 default and method=2 cc_fixed)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def _render_fig6_canary(canary_data: Mapping[str, Any], out_path: Path) -> None:
    """Render optional Figure 6: canary wall time (and t_stop) over run sequence."""
    plt = _get_plt()
    series_list = list(canary_data.get("series") or [])
    if not series_list:
        return
    n = len(series_list)
    fig, axes = plt.subplots(n, 1, figsize=(9.5, max(3.2, 3.0 * n)), squeeze=False)
    for idx, ser in enumerate(series_list):
        ax = axes[idx][0]
        pts = ser["points"]
        xs = [int(p["index"]) for p in pts]
        ws = [p["wall_sec"] for p in pts]
        ts = [p["t_stop"] for p in pts]
        ax.plot(xs, [w if w is not None else np.nan for w in ws], "o-", color="C0", lw=1.5, label="wall_sec")
        if any(t is not None for t in ts):
            ax.plot(
                xs,
                [t if t is not None else np.nan for t in ts],
                "s--",
                mfc="none",
                color="C1",
                lw=1.2,
                label="t_stop",
            )
        med = ser.get("median_wall_s")
        if med is not None and med > 0:
            ax.axhline(float(med), color="0.3", ls="--", lw=1.0, label=f"median ({_fmt_time(med)} s)")
            ax.axhspan(float(med) * 0.9, float(med) * 1.1, color="0.75", alpha=0.25, label="+-10% band")
        if len(xs) <= 16:
            ax.set_xticks(xs)
        ax.set_xlabel("canary run index", fontsize=8.5)
        ax.set_ylabel("time (s)", fontsize=8.5)
        ax.set_title(f"{ser['entry']} on {ser['instance']}", fontsize=9.5)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7.5, loc="best")
    fig.suptitle("Figure 6: Canary Drift Check", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def render_figures(result: Mapping[str, Any], out_dir: str | Path) -> dict[str, str]:
    """
    Render Figures 1-5 (and optional Figure 6 when canary data exists) into `out_dir`.

    Reads only `result["figures"]` (no access to raw records or manifest).

    Args:
        result: Dict returned by `analyze(...)`.
        out_dir: Output directory where PNG files are written.

    Returns:
        Mapping from figure key (`fig1`..`fig5`, optionally `fig6`) to PNG filename relative to `out_dir`.
    """
    odir = Path(out_dir)
    odir.mkdir(parents=True, exist_ok=True)
    figs = result.get("figures") or {}
    color_map = _entry_color_map(list(result.get("entries") or []))
    written: dict[str, str] = {}

    _render_fig1(figs.get("fig1") or {}, odir / "fig1_cactus.png", color_map=color_map)
    written["fig1"] = "fig1_cactus.png"
    _render_fig2(figs.get("fig2") or {}, odir / "fig2_scatter.png")
    written["fig2"] = "fig2_scatter.png"
    _render_fig3(figs.get("fig3") or {}, odir / "fig3_t1_t2_times.png", color_map=color_map)
    written["fig3"] = "fig3_t1_t2_times.png"
    _render_fig4(figs.get("fig4") or {}, odir / "fig4_t3_brackets.png", color_map=color_map)
    written["fig4"] = "fig4_t3_brackets.png"
    _render_fig5(figs.get("fig5") or {}, odir / "fig5_scaling.png")
    written["fig5"] = "fig5_scaling.png"

    fig6_data = figs.get("fig6") or {}
    if fig6_data.get("has_data"):
        _render_fig6_canary(fig6_data, odir / "fig6_canary.png")
        written["fig6"] = "fig6_canary.png"

    return written


# ======================================================================================================================
# Markdown report renderer (all lines <= 120 columns, no hard-coded verdicts)
# ======================================================================================================================

def _wrap(text: str, width: int = 120, initial_indent: str = "", subsequent_indent: str = "") -> list[str]:
    """Wrap prose to at most `width` columns cleanly."""
    if not text:
        return [""]
    return textwrap.wrap(
        text,
        width=width,
        initial_indent=initial_indent,
        subsequent_indent=subsequent_indent,
        break_long_words=True,
        break_on_hyphens=False,
    )


def _format_md_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    aligns: Sequence[str],
) -> list[str]:
    """Render a compact Markdown table without excess padding."""
    align_map = {"l": ":---", "r": "---:", "c": ":---:"}
    out = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(align_map.get(a, ":---") for a in aligns) + " |",
    ]
    for row in rows:
        safe = [str(c).replace("\n", " ").replace("|", "/") for c in row]
        out.append("| " + " | ".join(safe) + " |")
    return out


def _max_table_width(lines: Sequence[str]) -> int:
    """Return the maximum line length in a rendered Markdown table."""
    return max((len(ln) for ln in lines), default=0)


def _truncate_cell(text: str, max_len: int) -> str:
    """Truncate a Markdown table cell cleanly to `max_len` characters."""
    s = str(text).replace("\n", " ").replace("|", "/")
    if len(s) <= max_len:
        return s
    if max_len <= 3:
        return s[:max_len]
    return s[: max_len - 3] + "..."


def _render_table_a_md(
    result: Mapping[str, Any],
    entry_codes: Mapping[str, str],
) -> list[str]:
    """Render Table A (per tier and all tiers combined) in Markdown <= 120 cols."""
    out: list[str] = ["## Table A: Summary by Tier and Engine", ""]
    table_a = result.get("table_a") or {}
    tiers = list((result.get("meta") or {}).get("tiers") or TIER_ORDER_DEFAULT)
    scopes = tiers + ["all"]

    for scope in scopes:
        tinfo = table_a.get(scope) or {}
        budget = tinfo.get("budget_s")
        man_b = tinfo.get("budget_manifest_s")
        if scope == "all":
            heading = f"### All tiers ({tinfo.get('instances_manifest', 0)} manifest instances)"
        else:
            if budget is not None and man_b is not None and abs(float(budget) - float(man_b)) > 1e-6:
                b_str = f", budget {_fmt_time(budget)} s (manifest {_fmt_time(man_b)} s)"
            elif budget is not None:
                b_str = f", budget {_fmt_time(budget)} s"
            else:
                b_str = ""
            heading = f"### {scope} ({tinfo.get('instances_manifest', 0)} manifest instances{b_str})"
        out += [heading, ""]
        if not tinfo.get("has_data"):
            out += ["no data", ""]
            continue

        headers = [
            "Code",
            "Cert/Tot",
            "Inst",
            "PAR-2 (s)",
            "P2 ratio",
            "All-sol med",
            "Speedup [95% CI]",
            "W/T/L (+only)",
            "Marg. VBS",
        ]
        aligns = ["l", "r", "r", "r", "r", "r", "r", "r", "r"]
        all_rows = list(tinfo.get("rows") or [])
        if tinfo.get("vbs_row"):
            all_rows.append(tinfo["vbs_row"])

        md_rows: list[list[str]] = []
        missing_notes: list[str] = []
        for r in all_rows:
            entry = str(r["entry"])
            code = "VBS" if r.get("is_vbs") else entry_codes.get(entry, entry)
            if r.get("is_baseline"):
                code = f"{code}*"
            ct = f"{r['runs_certified']}/{r['runs_total']}"
            inst_s = f"{r['instances_certified']}/{r['instances_with_data']}({r['instances_manifest']})"
            if int(r.get("instances_missing") or 0) > 0 and not r.get("is_vbs"):
                missing_notes.append(f"{code}: {r['instances_missing']} missing")
            p2 = f"{_fmt_time(r['par2_mean'])}({r['par2_n']})" if r["par2_mean"] is not None else "no data"
            p2r = _fmt_ratio(r["par2_ratio_vs_baseline"])
            cmed = (
                f"{_fmt_time(r['common_all_solved_median_s'])}({r['common_all_solved_n']})"
                if r["common_all_solved_median_s"] is not None
                else "-"
            )
            if r["speedup_geomean"] is not None:
                sp_s = (
                    f"{_fmt_ratio(r['speedup_geomean'])} "
                    f"[{_fmt_ratio(r['speedup_ci_low'])},{_fmt_ratio(r['speedup_ci_high'])}]"
                    f"({r['speedup_n']})"
                )
            else:
                sp_s = "-"
            if r["wins"] is not None:
                wtl_s = f"{r['wins']}/{r['ties']}/{r['losses']}"
                oe = int(r.get("only_entry") or 0)
                ob = int(r.get("only_baseline") or 0)
                if oe or ob:
                    wtl_s += f"(+{oe}/-{ob})"
            else:
                wtl_s = "-"
            mv = f"+{_fmt_time(r['marginal_vbs'])}" if r.get("marginal_vbs") is not None else "-"
            md_rows.append([code, ct, inst_s, p2, p2r, cmed, sp_s, wtl_s, mv])

        tlines = _format_md_table(headers, md_rows, aligns)
        # Ensure every Table A row stays <= 120 chars even with huge numbers.
        if _max_table_width(tlines) > 120:
            compact_rows = []
            for row in md_rows:
                r2 = list(row)
                r2[6] = r2[6].replace(" ", "")
                compact_rows.append(r2)
            tlines = _format_md_table(headers, compact_rows, aligns)
        out.extend(tlines)
        out.append("")
        if missing_notes:
            out.extend(_wrap("Missing instances (not scored in PAR-2): " + "; ".join(missing_notes) + "."))
            out.append("")

    out.extend(
        _wrap(
            "Notes: `*` marks the baseline entry. `Inst` is `solved-by-median / with-data (manifest)`. "
            "`P2 ratio` is `PAR-2(baseline) / PAR-2(entry)` on common instances (>1 favors entry). "
            "`All-sol med` is the median of per-instance median times on instances solved by all entries in "
            "the tier. `Speedup` is the geometric mean of `t_base / t_entry` on both-solved instances with "
            "bootstrap 95% CI. `W/T/L` counts wins/ties/losses on both-solved instances (tie: within +-10% "
            "or overlapping seed min-max ranges), plus `(+only_entry/-only_base)` when nonzero. "
            "`Marg. VBS` is `PAR-2(VBS without entry) - PAR-2(VBS)`."
        )
    )
    out.append("")
    return out


def _render_table_b_md(
    result: Mapping[str, Any],
    entry_codes: Mapping[str, str],
) -> list[str]:
    """Render Table B (per instance, grouped by tier) chunking wide entry sets to stay <= 120 cols."""
    out: list[str] = ["## Table B: Per-Instance Results and VBS", ""]
    table_b = result.get("table_b") or {}
    tiers = list((result.get("meta") or {}).get("tiers") or TIER_ORDER_DEFAULT)

    for tier in tiers:
        tinfo = table_b.get(tier) or {}
        budget = tinfo.get("budget_s")
        man_b = tinfo.get("budget_manifest_s")
        if budget is not None and man_b is not None and abs(float(budget) - float(man_b)) > 1e-6:
            b_str = f" (budget {_fmt_time(budget)} s; manifest {_fmt_time(man_b)} s)"
        elif budget is not None:
            b_str = f" (budget {_fmt_time(budget)} s)"
        else:
            b_str = ""
        out += [f"### {tier}{b_str}", ""]
        if not tinfo.get("has_data"):
            out += ["no data", ""]
            continue

        entries = list(tinfo.get("entries") or [])
        rows = list(tinfo.get("rows") or [])

        # Prepare display columns: each entry + VBS.
        col_defs: list[tuple[str, str]] = [(entry_codes.get(e, e), e) for e in entries]
        col_defs.append(("VBS", "__VBS__"))

        inst_labels = [_short_inst_label(str(r["instance"]), 46) for r in rows]
        ref_labels = [str(r["ref_display"]) for r in rows]

        # Build formatted cells per column key.
        col_vals: dict[str, list[str]] = {}
        for _, key in col_defs:
            vals: list[str] = []
            for r in rows:
                if key == "__VBS__":
                    v = r["vbs"]
                    if not v["has_data"]:
                        vals.append("no data")
                    elif v["solved_by_median"]:
                        bcodes = ",".join(entry_codes.get(be, be) for be in v["best_entries"][:2])
                        vals.append(f"{_fmt_time(v['med_par2'])} ({bcodes})")
                    elif v["best_lb"] is not None or v["best_ub"] is not None:
                        lb_s = "-" if v["best_lb"] is None else str(int(v["best_lb"]))
                        ub_s = "-" if v["best_ub"] is None else str(int(v["best_ub"]))
                        conflict_flag = (
                            "!"
                            if (v["best_lb"] is not None and v["best_ub"] is not None and v["best_lb"] > v["best_ub"])
                            else ""
                        )
                        vals.append(f"[{lb_s}, {ub_s}]{conflict_flag}")
                    else:
                        vals.append(_fmt_time(v["med_par2"]))
                else:
                    txt = str((r.get("formatted") or {}).get(key, "no data"))
                    if key in (r.get("best_entries") or []) and txt != "no data":
                        txt = f"**{txt}**"
                    vals.append(txt)
            col_vals[key] = vals

        # Greedily pack columns into chunks so every line of every Markdown table is <= 120 chars.
        chunks: list[list[tuple[str, str]]] = []
        cur_chunk: list[tuple[str, str]] = []

        def _chunk_fits(chunk: Sequence[tuple[str, str]], include_ref: bool) -> bool:
            hdrs = ["Instance"] + (["Ref"] if include_ref else []) + [h for h, _ in chunk]
            md_r = []
            for idx in range(len(rows)):
                rr = [inst_labels[idx]] + ([ref_labels[idx]] if include_ref else [])
                rr.extend(col_vals[k][idx] for _, k in chunk)
                md_r.append(rr)
            lines = _format_md_table(hdrs, md_r, ["l"] * len(hdrs))
            return _max_table_width(lines) <= 120

        for col in col_defs:
            include_ref = len(chunks) == 0
            cand = cur_chunk + [col]
            if cur_chunk and not _chunk_fits(cand, include_ref):
                chunks.append(cur_chunk)
                cur_chunk = [col]
            else:
                cur_chunk = cand
        if cur_chunk:
            chunks.append(cur_chunk)

        for cidx, chunk in enumerate(chunks):
            include_ref = cidx == 0
            if len(chunks) > 1:
                col_list = ", ".join(h for h, _ in chunk)
                out.append(f"*{tier} part {cidx + 1}/{len(chunks)} ({col_list}):*")
                out.append("")
            hdrs = ["Instance"] + (["Ref"] if include_ref else []) + [h for h, _ in chunk]
            aligns = ["l"] + (["c"] if include_ref else []) + ["r"] * len(chunk)
            md_r: list[list[str]] = []
            for idx in range(len(rows)):
                rr = [inst_labels[idx]] + ([ref_labels[idx]] if include_ref else [])
                rr.extend(col_vals[k][idx] for _, k in chunk)
                md_r.append(rr)
            out.extend(_format_md_table(hdrs, md_r, aligns))
            out.append("")

        # Marginal VBS contribution summary for this tier.
        marg = tinfo.get("marginal_vbs") or {}
        if marg:
            m_parts = [
                f"{entry_codes.get(e, e)}=+{_fmt_time(marg[e])} s"
                for e in entries
                if marg.get(e) is not None
            ]
            vbs_p2 = tinfo.get("vbs_par2")
            summary = (
                f"{tier} VBS PAR-2: {_fmt_time(vbs_p2)} s. "
                f"Marginal VBS contributions: {', '.join(m_parts)}."
            )
            out.extend(_wrap(summary))
            out.append("")

    return out


def _render_table_c_md(
    result: Mapping[str, Any],
    load_meta: Sequence[Mapping[str, Any]] | None = None,
) -> list[str]:
    """Render Table C (integrity checks) in Markdown <= 120 cols; prints 'none' when empty."""
    tc = result.get("table_c") or {}
    out: list[str] = ["## Table C: Integrity Audit", ""]

    # 1. Wrong results
    out += ["### Wrong results (`status == wrong`)", ""]
    wrong = tc.get("wrong") or []
    if not wrong:
        out += ["none", ""]
    else:
        rows = []
        notes: list[str] = []
        for r in wrong:
            rid = str(r.get("run_id") or "")[:10]
            vtext = "; ".join(r.get("violations") or []) or "-"
            vshort = _truncate_cell(vtext, 28)
            if vshort != vtext:
                notes.append(f"`{rid}`: {vtext}")
            rows.append(
                [
                    rid,
                    _truncate_cell(str(r.get("entry") or ""), 24),
                    _truncate_cell(str(r.get("instance") or ""), 28),
                    f"[{r.get('lb')},{r.get('ub')}]",
                    vshort,
                ]
            )
        out.extend(_format_md_table(["run_id", "entry", "instance", "[lb,ub]", "violations"], rows, ["l"] * 5))
        out.append("")
        for note in notes:
            out.extend(_wrap(note, initial_indent="- ", subsequent_indent="  "))
        if notes:
            out.append("")

    # 2. Cross-run LB conflicts
    out += ["### Cross-run lower-bound conflicts (`lb > min verified witness weight`)", ""]
    lbc = tc.get("lb_conflicts") or []
    if not lbc:
        out += ["none", ""]
    else:
        rows = []
        for r in lbc:
            rows.append(
                [
                    str(r.get("run_id") or "")[:10],
                    _truncate_cell(str(r.get("entry") or ""), 22),
                    _truncate_cell(str(r.get("instance") or ""), 30),
                    str(r.get("lb")),
                    str(r.get("min_verified_weight")),
                    str(r.get("witness_run_id") or "")[:10],
                ]
            )
        out.extend(
            _format_md_table(
                ["run_id", "entry", "instance", "lb", "min_w", "witness_run"],
                rows,
                ["l", "l", "l", "r", "r", "l"],
            )
        )
        out.append("")

    # 3. Certified runs with unverified upper bound
    out += ["### Certified runs with unverified upper bound (`ub_verified != true`)", ""]
    unv = tc.get("unverified_certified") or []
    if not unv:
        out += ["none", ""]
    else:
        rows = []
        notes = []
        for r in unv:
            rid = str(r.get("run_id") or "")[:10]
            reason = str(r.get("unverified_reason") or "not verified")
            rshort = _truncate_cell(reason, 24)
            if rshort != reason:
                notes.append(f"`{rid}`: {reason}")
            rows.append(
                [
                    rid,
                    _truncate_cell(str(r.get("entry") or ""), 22),
                    _truncate_cell(str(r.get("instance") or ""), 28),
                    str(r.get("ub")),
                    str(r.get("ub_verified")),
                    rshort,
                ]
            )
        out.extend(
            _format_md_table(
                ["run_id", "entry", "instance", "ub", "verified", "reason"],
                rows,
                ["l", "l", "l", "r", "c", "l"],
            )
        )
        out.append("")
        for note in notes:
            out.extend(_wrap(note, initial_indent="- ", subsequent_indent="  "))
        if notes:
            out.append("")

    # 4. Contaminated runs
    excl_str = "excluded from main tables/figs" if tc.get("exclude_contaminated") else "included in main tables/figs"
    out += [f"### Contaminated runs (`env.contaminated == true`, {excl_str})", ""]
    contam = tc.get("contaminated") or []
    if not contam:
        out += ["none", ""]
    else:
        rows = []
        for r in contam:
            rows.append(
                [
                    str(r.get("run_id") or "")[:10],
                    _truncate_cell(str(r.get("entry") or ""), 22),
                    _truncate_cell(str(r.get("instance") or ""), 28),
                    str(r.get("threads") or "-"),
                    _fmt_time(r.get("wall_sec")),
                    _fmt_time(r.get("ext_cpu_s")),
                    _fmt_time(r.get("steal_s")),
                ]
            )
        out.extend(
            _format_md_table(
                ["run_id", "entry", "instance", "p", "wall(s)", "ext_cpu", "steal"],
                rows,
                ["l", "l", "l", "r", "r", "r", "r"],
            )
        )
        out.append("")

    # 5. Errors and excluded/bad records
    out += ["### Error runs and excluded records", ""]
    errs = tc.get("errors") or []
    excl = tc.get("excluded") or []
    bad_load = [m for m in (load_meta or []) if int(m.get("bad_lines") or 0) > 0]
    if not errs and not excl and not bad_load:
        out += ["none", ""]
    else:
        if errs:
            rows = []
            notes = []
            for r in errs:
                rid = str(r.get("run_id") or "")[:10]
                msg = str(r.get("msg") or "-")
                mshort = _truncate_cell(msg, 32)
                if mshort != msg:
                    notes.append(f"`{rid}`: {msg}")
                rows.append(
                    [
                        rid,
                        _truncate_cell(str(r.get("entry") or ""), 22),
                        _truncate_cell(str(r.get("instance") or ""), 26),
                        str(r.get("kind") or "-"),
                        mshort,
                    ]
                )
            out.extend(_format_md_table(["run_id", "entry", "instance", "kind", "msg"], rows, ["l"] * 5))
            out.append("")
            for note in notes:
                out.extend(_wrap(note, initial_indent="- ", subsequent_indent="  "))
            if notes:
                out.append("")
        if excl:
            out.append("Excluded records (manifest mismatch):")
            for r in excl:
                txt = f"`{str(r.get('run_id') or '')[:10]}` ({r.get('entry')}, {r.get('instance')}): {r.get('reason')}"
                out.extend(_wrap(txt, initial_indent="- ", subsequent_indent="  "))
            out.append("")
        if bad_load:
            out.append("Skipped malformed JSON lines during session load:")
            for m in bad_load:
                txt = f"{m.get('session_id')}: {m.get('bad_lines')} bad line(s)"
                out.extend(_wrap(txt, initial_indent="- ", subsequent_indent="  "))
            out.append("")
    return out


def _kv_table_rows(pairs: Sequence[tuple[str, str]], max_val_len: int = 92) -> list[list[str]]:
    """Build key-value table rows wrapping long values across continuation rows to stay <= 120 cols."""
    rows: list[list[str]] = []
    for key, val in pairs:
        v = str(val if val is not None else "-")
        chunks = _wrap(v, width=max_val_len) if len(v) > max_val_len else [v]
        for idx, ch in enumerate(chunks):
            rows.append([key if idx == 0 else "", ch])
    return rows


def _render_table_d_md(result: Mapping[str, Any]) -> list[str]:
    """Render Table D (environment from env.json and engine versions from records) <= 120 cols."""
    td = result.get("table_d") or {}
    out: list[str] = ["## Table D: Environment and Engine Versions", ""]

    envs = td.get("envs") or []
    missing_eids = td.get("missing_env_ids") or []
    out += ["### Hardware and software environments (`env.json`)", ""]
    if not envs:
        out += ["no data", ""]
    else:
        for env in envs:
            eid = str(env.get("env_id") or "")
            out += [f"#### Environment `{eid[:12]}` ({env.get('host') or 'unknown'})", ""]
            cpu = env.get("cpu") or {}
            topo = (
                f"{cpu.get('sockets')}x{cpu.get('cores_per_socket')}x{cpu.get('threads_per_core')} "
                f"({cpu.get('cpus')} logical CPUs)"
                if cpu.get("cpus") is not None
                else "-"
            )
            mem_gib = env.get("memory_gib")
            mem_s = f"{mem_gib:.1f} GiB" if mem_gib is not None else "-"
            comp = env.get("compiler") or {}
            py = env.get("python") or {}
            qc = env.get("qec_cert") or {}
            dm = env.get("dist_m4ri") or {}
            harn = env.get("harness") or {}
            man = env.get("manifest") or {}
            inp = env.get("inputs") or {}

            qc_s = (
                f"{qc.get('version')} commit={str(qc.get('commit') or '-')[:10]} "
                f"dirty={qc.get('dirty')} lib={str(qc.get('lib_sha256') or '-')[:12]}"
                if qc.get("version") or qc.get("commit")
                else "not recorded"
            )
            dm_s = (
                f"bin={dm.get('version')} commit={str(dm.get('commit') or '-')[:10]} "
                f"dirty={dm.get('dirty')} sha={str(dm.get('sha256') or '-')[:12]}; "
                f"py={dm.get('py_version')} sha={str(dm.get('py_sha256') or '-')[:12]}"
                if dm.get("version") or dm.get("commit")
                else "not recorded"
            )
            py_s = (
                f"Python {py.get('version') or '-'}, stim={py.get('stim') or '-'}, "
                f"numpy={py.get('numpy') or '-'}, scipy={py.get('scipy') or '-'}, "
                f"matplotlib={py.get('matplotlib') or '-'}"
            )
            inp_s = (
                f"checked={inp.get('checked')}, files={inp.get('files')}, "
                f"sha256_ok={inp.get('sha256_ok')}, warmed={inp.get('page_cache_warmed')}"
                if any(v is not None for v in inp.values())
                else "-"
            )
            numa_raw = cpu.get("numa")
            if isinstance(numa_raw, Mapping):
                numa_s = "; ".join(f"{k}={v}" for k, v in sorted(numa_raw.items())) or "-"
            else:
                numa_s = str(numa_raw) if numa_raw is not None else "-"
            flags_raw = comp.get("flags")
            if isinstance(flags_raw, Mapping):
                flags_s = ", ".join(f"{k}={v}" for k, v in sorted(flags_raw.items())) or "-"
            else:
                flags_s = str(flags_raw) if flags_raw is not None else "-"
            pairs = [
                ("env_id", eid),
                ("created / host", f"{env.get('created') or '-'} / {env.get('host') or '-'}"),
                ("OS / kernel", f"{env.get('os') or '-'} / {env.get('kernel') or '-'}"),
                ("VM / sched", f"{env.get('vm') or '-'} / autogroup={env.get('sched_autogroup')}"),
                ("CPU model", str(cpu.get("model") or "-")),
                ("CPU topology", f"{topo}; numa={numa_s}; flags={str(cpu.get('flags_sha1') or '-')[:12]}"),
                ("Memory", mem_s),
                ("Compiler", f"{comp.get('cc') or '-'}; flags={flags_s}"),
                ("libm4ri", str(env.get("m4ri") or "-")),
                ("Python libs", py_s),
                ("qec_cert", qc_s),
                ("dist_m4ri", dm_s),
                ("Harness", str(harn.get("version") or "-")),
                (
                    "Manifest",
                    f"{man.get('id') or '-'}/{man.get('version') or '-'} sha={str(man.get('sha256') or '-')[:12]}",
                ),
                ("Inputs", inp_s),
            ]
            out.extend(_format_md_table(["Field", "Value"], _kv_table_rows(pairs, 92), ["l", "l"]))
            out.append("")

    if missing_eids:
        out.extend(
            _wrap(
                "Note: records reference env_id(s) without a loaded env.json: "
                + ", ".join(f"`{e[:12]}`" for e in missing_eids)
                + "."
            )
        )
        out.append("")

    # Engine versions from run records.
    out += ["### Engine versions (from run records)", ""]
    engs = td.get("engines") or []
    if not engs:
        out += ["no data", ""]
    else:
        rows = []
        for e in engs:
            flag = " [MULTI-VER]" if e.get("version_conflict") else ""
            ver_s = str(e.get("version") or "-")
            com_s = str(e.get("commit") or "-")[:9] + ("*" if e.get("dirty") else "")
            sha_s = str(e.get("sha256") or "-")[:12]
            deps = e.get("deps") or {}
            dep_s = ", ".join(f"{k}={v}" for k, v in sorted(deps.items())) if deps else "-"
            rows.append(
                [
                    _truncate_cell(f"{e['entry']}{flag}", 34),
                    _truncate_cell(ver_s, 10),
                    com_s,
                    sha_s,
                    str(e.get("n_runs") or 0),
                    _truncate_cell(dep_s, 26),
                ]
            )
        out.extend(
            _format_md_table(
                ["Entry", "Version", "Commit", "sha256[:12]", "Runs", "Deps"],
                rows,
                ["l", "l", "l", "l", "r", "l"],
            )
        )
        out.append("")
        if td.get("multi_version_entries"):
            out.extend(
                _wrap(
                    "WARNING: multiple engine versions/builds detected within the same entry: "
                    + ", ".join(td["multi_version_entries"])
                    + "."
                )
            )
            out.append("")
    return out


def _render_canary_md(result: Mapping[str, Any], fig_files: Mapping[str, str]) -> list[str]:
    """Render optional Canary Drift section."""
    canary = result.get("canary") or {}
    out: list[str] = ["## Canary Drift Check", ""]
    if not canary.get("has_data"):
        out += ["no data", ""]
        return out
    rows = []
    for s in canary.get("series") or []:
        rows.append(
            [
                _truncate_cell(str(s["entry"]), 28),
                _truncate_cell(str(s["instance"]), 30),
                f"{s['n_certified']}/{s['n']}",
                _fmt_time(s["median_wall_s"]),
                _fmt_time(s["min_wall_s"]),
                _fmt_time(s["max_wall_s"]),
                _fmt_num(s["spread_pct"], 1) + "%" if s["spread_pct"] is not None else "-",
            ]
        )
    out.extend(
        _format_md_table(
            ["Entry", "Instance", "Cert/N", "Med(s)", "Min(s)", "Max(s)", "Spread"],
            rows,
            ["l", "l", "r", "r", "r", "r", "r"],
        )
    )
    out.append("")
    if "fig6" in fig_files:
        out += [f"![Figure 6: Canary Drift]({fig_files['fig6']})", ""]
        out.extend(_wrap("*" + str((result.get("captions") or {}).get("fig6") or "") + "*"))
        out.append("")
    return out


def _render_scaling_table_md(result: Mapping[str, Any]) -> list[str]:
    """Render a compact summary table for Figure 5 scaling series."""
    fig5 = (result.get("figures") or {}).get("fig5") or {}
    if not fig5.get("has_data"):
        return ["no data", ""]
    rows = []
    for ser in fig5.get("series") or []:
        iid = _short_inst_label(str(ser["instance"]), 32)
        cfg = str(ser["config"])
        for pt in ser.get("points") or []:
            sp = _fmt_num(pt["speedup"], 2) if pt.get("speedup") is not None else "-"
            eff = _fmt_num(pt["efficiency"], 2) if pt.get("efficiency") is not None else "-"
            rows.append(
                [
                    iid,
                    cfg,
                    str(pt["threads"]),
                    str(pt["n"]),
                    _fmt_time(pt["med_wall_s"]),
                    f"[{_fmt_time(pt['min_wall_s'])}, {_fmt_time(pt['max_wall_s'])}]",
                    sp,
                    eff,
                ]
            )
    return _format_md_table(
        ["Instance", "Config", "p", "N", "Med(s)", "[Min, Max]", "Speedup", "Eff"],
        rows,
        ["l", "l", "r", "r", "r", "r", "r", "r"],
    ) + [""]


def render_report(
    result: Mapping[str, Any],
    fig_files: Mapping[str, str] | None = None,
    *,
    load_meta: Sequence[Mapping[str, Any]] | None = None,
) -> str:
    """
    Render the full Markdown report (`REPORT.md`) from an `analyze(...)` result dict.

    Guarantees every line is at most 120 columns and every statement is computed from `result`
    with zero hard-coded verdicts.

    Args:
        result: Dict returned by `analyze(...)`.
        fig_files: Optional mapping `{fig1..fig6: relative_png_filename}`.
        load_meta: Optional session loader metadata list (for reporting skipped malformed lines).

    Returns:
        Markdown string ending with a newline.
    """
    figs = dict(fig_files or {})
    meta = result.get("meta") or {}
    caps = result.get("captions") or {}
    entry_codes = result.get("entry_codes") or {}
    entries = list(result.get("entries") or [])
    baseline = str(meta.get("baseline") or DEFAULT_BASELINE)

    out: list[str] = [
        f"# Certification Benchmark Report ({meta.get('manifest_id')}/{meta.get('manifest_version')})",
        "",
    ]

    # Computed provenance summary (no hard-coded verdicts).
    summary_line = (
        f"Manifest: `{meta.get('manifest_id')}` v{meta.get('manifest_version')}. "
        f"Baseline: `{baseline}`. "
        f"Run dates: {meta.get('date_range') or 'unknown'}. "
        f"Records: {meta.get('records_valid', 0)} deduplicated "
        f"({meta.get('records_main', 0)} main default-league, {meta.get('records_canary', 0)} canary, "
        f"{meta.get('records_duplicates_overwritten', 0)} duplicate run_id overwritten, "
        f"{meta.get('records_excluded', 0)} excluded). "
        f"Contaminated: {meta.get('contaminated_total', 0)} total "
        f"({'excluded' if meta.get('exclude_contaminated') else 'included'} in main tables/figures). "
        f"Bootstrap: {meta.get('n_boot', 0)} resamples (seed={meta.get('boot_seed', 1)})."
    )
    out.extend(_wrap(summary_line))
    out.append("")
    out.extend(_wrap(WALL_TIME_SENTENCE))
    out.append("")

    # Entry key (maps E0, E1, ... to full engine entry names).
    out += ["## Engine Entries", ""]
    if not entries:
        out += ["no data", ""]
    else:
        eng_lookup: dict[str, dict[str, Any]] = {}
        for er in (result.get("table_d") or {}).get("engines") or []:
            eng_lookup.setdefault(str(er["entry"]), er)
        key_rows = []
        for e in entries:
            er = eng_lookup.get(e) or {}
            ver = str(er.get("version") or "-")
            com = str(er.get("commit") or "-")[:9] + ("*" if er.get("dirty") else "")
            role = "baseline" if e == baseline else ""
            key_rows.append([entry_codes.get(e, e), _truncate_cell(e, 38), ver, com, role])
        out.extend(
            _format_md_table(
                ["Code", "Entry", "Version", "Commit", "Role"],
                key_rows,
                ["l", "l", "l", "l", "l"],
            )
        )
        out.append("")

    # Table A
    out.extend(_render_table_a_md(result, entry_codes))

    # Figures 1-4
    out += ["## Certification Figures (Figures 1-4)", ""]
    fig_titles = [
        ("fig1", "Figure 1: Cactus / ECDF by Tier"),
        ("fig2", "Figure 2: Baseline vs. Entry Scatter"),
        ("fig3", "Figure 3: Per-Instance Time on T1-T2"),
        ("fig4", "Figure 4: T3 Brackets"),
    ]
    for fkey, ftitle in fig_titles:
        out += [f"### {ftitle}", ""]
        fdata = (result.get("figures") or {}).get(fkey) or {}
        if not fdata.get("has_data"):
            out += ["no data", ""]
        if fkey in figs:
            out += [f"![{ftitle}]({figs[fkey]})", ""]
        cap = str(caps.get(fkey) or "")
        if cap:
            out.extend(_wrap(f"*{cap}*"))
            out.append("")

    # Table B
    out.extend(_render_table_b_md(result, entry_codes))

    # Figure 5 (Thread scaling)
    out += ["## Figure 5: Thread Scaling (`dist_m4ri` binary)", ""]
    if not ((result.get("figures") or {}).get("fig5") or {}).get("has_data"):
        out += ["no data", ""]
    if "fig5" in figs:
        out += [f"![Figure 5: Thread Scaling]({figs['fig5']})", ""]
    cap5 = str(caps.get("fig5") or "")
    if cap5:
        out.extend(_wrap(f"*{cap5}*"))
        out.append("")
    if ((result.get("figures") or {}).get("fig5") or {}).get("has_data"):
        out.extend(_render_scaling_table_md(result))

    # Table C
    out.extend(_render_table_c_md(result, load_meta=load_meta))

    # Canary drift
    out.extend(_render_canary_md(result, figs))

    # Table D
    out.extend(_render_table_d_md(result))

    return "\n".join(out).rstrip() + "\n"


# ======================================================================================================================
# CLI
# ======================================================================================================================

def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for `tools/bench/analyze.py`."""
    parser = argparse.ArgumentParser(
        description=(
            "Analyze dist-m4ri certification benchmark sessions (dm4ri-bench-run/1) and write "
            "REPORT.md, summary.json, and fig1..fig5 PNGs."
        )
    )
    parser.add_argument("--manifest", required=True, help="Path to manifest JSON (e.g. certify_v1.json).")
    parser.add_argument(
        "--session",
        action="append",
        required=True,
        dest="sessions",
        help="Session directory (repeatable; later sessions win on duplicate run_id).",
    )
    parser.add_argument("--out", required=True, help="Output directory for REPORT.md, summary.json, and figures.")
    parser.add_argument(
        "--exclude-contaminated",
        action="store_true",
        help="Exclude runs with env.contaminated == true from main tables and figures (still listed in Table C).",
    )
    parser.add_argument(
        "--baseline",
        default=DEFAULT_BASELINE,
        help=f"Baseline engine entry (default: {DEFAULT_BASELINE}).",
    )
    parser.add_argument(
        "--n-boot",
        type=int,
        default=2000,
        help="Bootstrap resamples for geometric-mean speedup 95%% CI (default: 2000).",
    )
    parser.add_argument(
        "--boot-seed",
        type=int,
        default=1,
        help="Fixed RNG seed for bootstrap resampling (default: 1).",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip rendering PNG figures (still computes figure data and writes summary.json and REPORT.md).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point: load manifest and sessions, run `analyze`, render figures and REPORT.md."""
    args = _parse_args(argv)
    manifest_path = Path(args.manifest)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    records, envs, load_meta = load_sessions(args.sessions)
    result = analyze(
        manifest,
        records,
        envs,
        baseline=args.baseline,
        exclude_contaminated=bool(args.exclude_contaminated),
        n_boot=int(args.n_boot),
        boot_seed=int(args.boot_seed),
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig_files: dict[str, str] = {}
    if not args.no_figures:
        fig_files = render_figures(result, out_dir)

    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    report_md = render_report(result, fig_files=fig_files, load_meta=load_meta)
    report_path = out_dir / "REPORT.md"
    report_path.write_text(report_md, encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
