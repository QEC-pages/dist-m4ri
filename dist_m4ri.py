"""
dist_m4ri.py: Python wrapper for the multithreaded dist_m4ri distance calculator.

Provides high-level APIs for computing:
- Classical code distance: compute_classical_distance(...)
- Single-sided quantum code distance: compute_quantum_distance(...)
- CSS quantum code distance: compute_css_distance(...)
- Detector Error Model (DEM) distance: compute_dem_distance(...)

Supports distance caching, codeword export, and fallback/option for the codedistance library.
Since dist_m4ri natively handles multithreading (via POSIX threads and dynamic bracketing),
no Python-level threading or subprocess-per-core logic is necessary.
"""

import os
import sys
import json
import time
import random
import shutil
import hashlib
import tempfile
import threading
import subprocess
from pathlib import Path
from typing import List, Tuple, Union, Optional, Dict, Any, Set, Sequence, Callable

_codedistance_mod = None
_stim_mod = None

__version__ = "0.11.0"

# Debug bitmap (debug=N): bits 1..32768 are passed to the dist_m4ri binary (1: summary, 2: progress, 4: periodic
# status, 8: thread allocation and timing, 16: code parameters, 32: codeword supports, 64: arguments, 128: matrices,
# 256: codeword dump, 2048: no size cutoff; see its --morehelp), the higher bits are used by the Python wrapper only.
PY_DBG_COMMANDS = 1 << 16    # echo the dist_m4ri command lines and their run times
PY_DBG_CACHE = 1 << 17       # cache messages
PY_DBG_CIRCUITS = 1 << 18    # Stim circuit to DEM conversion details
PY_DBG_KEEP_FILES = 1 << 19  # keep (and list) the temporary matrix, codeword, and DEM files
BIN_DBG_MASK = 0xFFFF        # the debug bits passed to the dist_m4ri binary
VERBOSE_DEBUG = 3 | PY_DBG_COMMANDS | PY_DBG_CACHE  # debug bits implied by verbose=True (--verbose)


def parse_debug_value(val: Union[str, int]) -> int:
    """
    Parses a debug bitmap value: a non-negative decimal or hexadecimal (0x...) integer, as in the dist_m4ri binary
    (no octal: leading zeros are ignored).

    Args:
        val: The value, e.g., "3" or "0x10004", or an int.

    Returns:
        The debug bitmap.

    Raises:
        ValueError: If the value is not a non-negative integer.
    """
    msg = f"invalid debug={val!r}: expected a non-negative decimal or hexadecimal (0x...) integer"
    if isinstance(val, bool):
        raise ValueError(msg)
    if isinstance(val, int):
        value = val
    else:
        s = str(val).strip()
        try:
            value = int(s[2:], 16) if s.lower().startswith("0x") else int(s, 10)
        except ValueError:
            raise ValueError(msg) from None
    if value < 0:
        raise ValueError(msg)
    return value


def _get_codedistance():
    """Lazily imports the codedistance library only when requested."""
    global _codedistance_mod
    if _codedistance_mod is None:
        try:
            import codedistance
            _codedistance_mod = codedistance
        except ImportError:
            raise ImportError("codedistance library is requested but not installed.")
    return _codedistance_mod


def _get_stim():
    """Lazily imports stim only when requested."""
    global _stim_mod
    if _stim_mod is None:
        try:
            import stim
            _stim_mod = stim
        except ImportError:
            raise ImportError("stim library is requested but not installed.")
    return _stim_mod


def __getattr__(name: str) -> Any:
    """Lazily resolves module attributes without loading heavy dependencies at startup."""
    if name == "_HAS_STIM":
        try:
            import stim
            return True
        except ImportError:
            return False
    if name == "_HAS_CODEDISTANCE":
        try:
            import codedistance
            return True
        except ImportError:
            return False
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

# Global cache for distance results
_distance_cache: Dict[str, Any] = {}
_use_distance_cache: bool = True
_distance_cache_file: Optional[str] = None
_last_run_stats: Dict[str, Any] = {}
_last_css_stats: Dict[str, Dict[str, Any]] = {"X": {}, "Z": {}}

# Aliases for backward compatibility with vecdec.py
_css_distance_cache = _distance_cache
_use_css_distance_cache = _use_distance_cache


def set_distance_cache_file(filepath: Optional[Union[str, Path]] = None) -> None:
    """
    Sets the default JSON file for persistent distance caching.
    If the file exists, its contents are loaded into memory.
    """
    global _distance_cache_file
    if filepath is not None:
        _distance_cache_file = str(Path(filepath).resolve())
        load_distance_cache(_distance_cache_file)
    else:
        _distance_cache_file = None


def load_distance_cache(filepath: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """
    Loads distance cache from a JSON file into memory.
    Inspects cache version silently; if an incompatible version is detected in the future,
    triggers a warning and bypasses/updates the cache.
    """
    global _distance_cache, _distance_cache_file
    target_file = str(Path(filepath).resolve()) if filepath is not None else _distance_cache_file
    if target_file and os.path.isfile(target_file):
        try:
            with open(target_file, "r") as f:
                data = json.load(f)
            if isinstance(data, dict):
                cache_ver = data.pop("__version__", None)
                # If cache is from a future incompatible version, skip loading
                if cache_ver is not None and _parse_version(cache_ver) > _parse_version(__version__):
                    sys.stderr.write(
                        f"# Warning: Cache file '{target_file}' has newer version {cache_ver} "
                        f"(current {__version__}); ignoring incompatible cache.\n"
                    )
                    return _distance_cache
                # Sanitize any legacy CSS cache entries where dmax was set from dmin when dmax_X == dmax_Z == 0
                for k, v in data.items():
                    if isinstance(v, dict) and k.startswith("css:"):
                        if v.get("dmax_X", 0) == 0 and v.get("dmax_Z", 0) == 0 and "dmax_X" in v:
                            v["dmax"] = 0
                _distance_cache.update(data)
        except Exception as e:
            sys.stderr.write(f"# Warning: Failed to load distance cache from {target_file}: {e}\n")
    return _distance_cache


def save_distance_cache(filepath: Optional[Union[str, Path]] = None) -> None:
    """
    Saves the in-memory distance cache to a JSON file.
    Silently writes "__version__": __version__ into the file.
    Uses atomic write via a temporary file to prevent corruption.
    """
    global _distance_cache, _distance_cache_file
    target_file = str(Path(filepath).resolve()) if filepath is not None else _distance_cache_file
    if not target_file:
        return

    parent_dir = os.path.dirname(os.path.abspath(target_file)) or "."
    os.makedirs(parent_dir, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(suffix=".tmp", prefix="dist_cache_", dir=parent_dir)
    try:
        cache_data = dict(_distance_cache)
        cache_data["__version__"] = __version__
        with open(fd, "w") as f:
            json.dump(cache_data, f, indent=2)
        os.replace(temp_path, target_file)
    except Exception as e:
        if os.path.exists(temp_path):
            try: os.remove(temp_path)
            except OSError: pass
        sys.stderr.write(f"# Warning: Failed to save distance cache to {target_file}: {e}\n")


def clear_distance_cache(cache_file: Optional[Union[str, Path]] = None, clear_file: bool = False) -> None:
    """
    Clears all cached distance calculations from memory, and optionally deletes the persistent JSON file.
    """
    global _distance_cache, _distance_cache_file
    _distance_cache.clear()
    target_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    if clear_file and target_file and os.path.isfile(target_file):
        try:
            os.remove(target_file)
        except OSError:
            pass


def enable_distance_cache() -> None:
    """Enables distance caching."""
    global _use_distance_cache
    _use_distance_cache = True


def disable_distance_cache() -> None:
    """Disables distance caching."""
    global _use_distance_cache
    _use_distance_cache = False


def get_distance_cache() -> Dict[str, Any]:
    """Returns the global distance cache dictionary."""
    global _distance_cache
    return _distance_cache


def format_bounds_list(dmin: int, dmax: int, num_rw: int) -> List[int]:
    """
    Returns [dmin, dmax, num_rw] according to:
    - [d, d, 0] if known exactly
    - [dmin, 0, 0] if there is no upper bound (dmax == 0)
    - [0, dmax, num_rw] if there is no lower bound (dmin <= 1)
    - [dmin, dmax, num_rw] otherwise
    """
    eff_dmin = dmin if dmin > 1 else 0
    eff_dmax = dmax if dmax > 0 else 0
    eff_rw = num_rw if num_rw > 0 else 0
    if eff_dmin > 0 and eff_dmin == eff_dmax:
        return [eff_dmin, eff_dmax, 0]
    elif eff_dmax == 0:
        return [eff_dmin, 0, 0]
    elif eff_dmin == 0:
        return [0, eff_dmax, eff_rw]
    else:
        return [eff_dmin, eff_dmax, eff_rw]


def format_bounds_str(bounds: List[int]) -> str:
    """Formats a bounds list [dmin, dmax, num_rw] as a string, stating '(exact)' if bounds coincide."""
    dmin, dmax, num_rw = bounds[0], bounds[1], bounds[2]
    if dmin > 0 and dmin == dmax:
        return f"{dmin} {dmax} {num_rw} (exact)"
    return f"{dmin} {dmax} {num_rw}"


def rw_infoset_find_prob(n: int, rank: int, w: int) -> float:
    """
    Information-set estimate of the probability that a RW step finds a given codeword of weight w.

    Model: uniform random information sets (as in the dist_m4ri binary, rw_infoset_find_prob()).  A RW step
    reduces H (rank r) with a random column order; a codeword is found if exactly one of its w positions is among
    the k = n - r non-pivot columns: P1(w) = w C(n-w, k-1) / C(n, k).  Codewords of weight w > r + 1 are never
    found (P1 = 0).

    Args:
        n: Block length (number of columns of H).
        rank: Rank r of H.
        w: Codeword weight.

    Returns:
        P1(w).
    """
    import math
    k = n - rank
    if n <= 0 or rank < 0 or k <= 0 or w < 1 or w > n or w - 1 > rank:
        return 0.0

    def log_binom(a: int, b: int) -> float:
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

    return min(1.0, math.exp(math.log(w) + log_binom(n - w, k - 1) - log_binom(n, k)))


def rw_infoset_estimate(
    n: int, rank: int, steps: int, d: int, dmin: int = 1, step_time: float = 0.0, target: float = 0.01,
    steps_total: int = 0
) -> Optional[Dict[str, Any]]:
    """
    Information-set estimate of the probability that RW missed a codeword lighter than the upper bound d found.

    Among the weights max(dmin, 1) <= w < d (not excluded by CC), the weight hardest to find (the smallest P1(w),
    see rw_infoset_find_prob()) is used: a single codeword of this weight is missed in `steps` RW steps with uniform
    random permutations with probability (1 - P1)^steps.  The number of RW steps for the target miss probability is
    given in total, with the fraction steps / steps_total of uniform steps as in the run (localized windows, used in
    about half of the steps for n >= 500, are not counted).  Same as print_rw_infoset_estimate() in the binary.

    Args:
        n: Block length (number of columns of H).
        rank: Rank of H.
        steps: Number of completed RW steps with uniform random permutations.
        d: Smallest weight found (the upper bound dmax).
        dmin: Certified lower bound (smaller weights are excluded).
        step_time: Wall time per RW step in seconds (0: unknown).
        target: Target miss probability (default: 0.01).
        steps_total: Number of all completed RW steps (0: the same as `steps`).

    Returns:
        Dictionary with the weight `w`, the probability `p_find` = P1(w) per step, the miss probability `p_miss`,
        the total number of RW steps `steps_target` for the target miss probability (None if the codeword cannot be
        found), and the corresponding wall time `time_target` (None if unknown); or None if there are no weights to
        check.
    """
    import math
    w_lo, w_hi = max(dmin, 1), d - 1
    if w_hi < w_lo or n <= 0:
        return None
    w_hard, p1 = w_hi, rw_infoset_find_prob(n, rank, w_hi)  # the weight hardest to find (smallest P1)
    for w in range(w_lo, w_hi):
        pw = rw_infoset_find_prob(n, rank, w)
        if pw < p1:
            w_hard, p1 = w, pw
    steps_target: Optional[int]
    if p1 <= 0.0:
        p_miss, steps_target = 1.0, None
    elif p1 >= 1.0:
        p_miss, steps_target = (0.0 if steps > 0 else 1.0), 1
    else:
        lq = math.log1p(-p1)
        p_miss = math.exp(steps * lq) if steps > 0 else 1.0
        steps_target = math.ceil(math.log(target) / lq)
    if steps_target is not None and steps > 0 and steps_total > steps:
        steps_target = math.ceil(steps_target * steps_total / steps)
    time_target = steps_target * step_time if (steps_target is not None and step_time > 0.0) else None
    return {"w": w_hard, "p_find": p1, "p_miss": p_miss, "steps_target": steps_target, "time_target": time_target}


def _format_duration(t: float) -> str:
    """Formats a duration in seconds (s, min, h, or days), as in the dist_m4ri binary."""
    if t < 120.0:
        return f"{t:.3g} s"
    if t < 7200.0:
        return f"{t / 60.0:.3g} min"
    if t < 172800.0:
        return f"{t / 3600.0:.3g} h"
    return f"{t / 86400.0:.3g} days"


def _parse_stderr_stats(stderr: str) -> Dict[str, Any]:
    """Parses codeword and hit statistics from dist_m4ri stderr output."""
    import re
    stats: Dict[str, Any] = {"extra_weights": [], "rw_converged": False, "hit_warnings": []}
    if not stderr:
        return stats
    hit_warning: Optional[List[str]] = None
    skip_estimate = False
    for line in stderr.splitlines():
        line_s = line.strip()
        if line_s.startswith("#   information-set estimate"):  # printed by the binary after a hit warning
            skip_estimate = True
        elif skip_estimate and line_s.startswith("#   "):
            continue
        else:
            skip_estimate = False
        if skip_estimate:
            if hit_warning is not None:
                stats["hit_warnings"].append(" ".join(hit_warning))
                hit_warning = None
            continue
        if hit_warning is not None:
            if line_s.startswith("#   "):  # continuation of a hit uniformity warning
                hit_warning.append(line_s[1:].strip())
                continue
            stats["hit_warnings"].append(" ".join(hit_warning))
            hit_warning = None
        if line_s.startswith("# Warning: non-uniform RW hit counts"):
            hit_warning = [line_s[len("# Warning: "):]]
            continue
        if line_s.startswith("# RW information sets:"):
            m_rw = re.search(
                r"n=(\d+), rank\(H\)=(\d+), steps=(\d+) \(uniform permutations: (\d+)\), ksub=(\d+), "
                r"([0-9.eE+-]+) s/step per thread, (\d+) RW threads", line_s
            )
            if m_rw:
                stats["rw_n"] = int(m_rw.group(1))
                stats["rw_rank"] = int(m_rw.group(2))
                stats["rw_steps"] = int(m_rw.group(3))
                stats["rw_steps_uniform"] = int(m_rw.group(4))
                stats["rw_ksub"] = int(m_rw.group(5))
                stats["rw_step_time"] = float(m_rw.group(6))
                stats["rw_threads"] = int(m_rw.group(7))
            continue
        if "RW convergence reached" in line_s:  # progress line (debug & 2) or the stop reason (debug & 1)
            stats["rw_converged"] = True
        if line_s.startswith("# codewords accumulated: total=0"):
            stats["total_cws"] = 0
        elif line_s.startswith("# codewords accumulated:"):
            m = re.search(
                r"total=(\d+),\s*min_w=(\d+):\s*cws=(\d+),\s*total_hits=(\d+),\s*"
                r"hits min=(\d+),\s*max=(\d+),\s*avg=([0-9.]+),\s*stdev=([0-9.]+)",
                line_s
            )
            if m:
                stats["total_cws"] = int(m.group(1))
                stats["min_w"] = int(m.group(2))
                stats["cws"] = int(m.group(3))
                stats["total_hits"] = int(m.group(4))
                stats["hits_min"] = int(m.group(5))
                stats["hits_max"] = int(m.group(6))
                stats["hits_avg"] = float(m.group(7))
                stats["hits_stdev"] = float(m.group(8))
            m_set = re.search(
                r"<n>=([0-9.]+) over (\d+) cws of w=(\d+)\.\.(\d+) \(min_hits=(\d+), cov_cws=(-?\d+)\)", line_s
            )
            if m_set:
                stats["avg_hits"] = float(m_set.group(1))
                stats["set_cws"] = int(m_set.group(2))
                stats["set_w_lo"] = int(m_set.group(3))
                stats["set_w_hi"] = int(m_set.group(4))
                stats["min_hits"] = int(m_set.group(5))
                stats["cov_cws"] = int(m_set.group(6))
        elif line_s.startswith("# codewords w="):
            m_w = re.search(
                r"w=(\d+):\s*cws=(\d+),\s*total_hits=(\d+),\s*"
                r"hits min=(\d+),\s*max=(\d+),\s*avg=([0-9.]+),\s*stdev=([0-9.]+)",
                line_s
            )
            if m_w:
                stats["extra_weights"].append({
                    "w": int(m_w.group(1)),
                    "cws": int(m_w.group(2)),
                    "total_hits": int(m_w.group(3)),
                    "hits_min": int(m_w.group(4)),
                    "hits_max": int(m_w.group(5)),
                    "hits_avg": float(m_w.group(6)),
                    "hits_stdev": float(m_w.group(7)),
                })
    if hit_warning is not None:
        stats["hit_warnings"].append(" ".join(hit_warning))
    return stats


def explain_bounds(
    bounds: List[int],
    method: Optional[int] = None,
    label: str = "",
    stats: Optional[Dict[str, Any]] = None
) -> str:
    """
    Returns a human-readable explanation of [dmin, dmax, num_rw] following README.md.

    Args:
        bounds: [dmin, dmax, rw_steps] list.
        method: Optional solver method (1=RW, 2=CC, 3=Bracketing).
        label: Optional prefix/label (e.g. "dX", "dZ", "").
        stats: Optional dictionary of codeword and hit statistics parsed from dist_m4ri.

    Returns:
        Multi-line formatted explanation string.
    """
    dmin, dmax, num_rw = bounds[0], bounds[1], bounds[2]
    lines = []
    prefix = f"{label} " if label else ""

    # Lower bound explanation
    if dmin > 0 and dmin == dmax:
        lines.append(f"  {prefix}Lower bound (dmin = {dmin}): Exact distance certified (dmin == dmax == {dmin}).")
    elif dmin > 1:
        lines.append(f"  {prefix}Lower bound (dmin = {dmin}): All cluster weights w <= {dmin - 1} "
                     f"were exhaustively analyzed by CC without finding any non-trivial codewords.")
    else:
        lines.append(f"  {prefix}Lower bound (dmin = {dmin}): No non-trivial lower bound certified (dmin <= 1).")

    # Upper bound explanation
    if dmax > 0:
        lines.append(f"  {prefix}Upper bound (dmax = {dmax}): Weight of the smallest non-trivial codeword discovered.")
    else:
        lines.append(f"  {prefix}Upper bound (dmax = {dmax}): No non-trivial codeword discovered yet (dmax = 0).")

    # RW steps explanation (and why it is zero if num_rw == 0)
    if num_rw > 0:
        lines.append(f"  {prefix}Random window steps (rw_steps = {num_rw}): {num_rw} completed random "
                     f"information set searches across worker threads.")
    else:
        if dmin > 0 and dmin == dmax:
            lines.append(f"  {prefix}Random window steps (rw_steps = 0): Set to 0 because the exact distance "
                     f"d = {dmin} was proven by Connected Cluster search or certified bounds coincided.")
        elif method == 2:
            lines.append(f"  {prefix}Random window steps (rw_steps = 0): Set to 0 because Method 2 (Connected Cluster) "
                     f"is an exhaustive search that does not perform random information set (RW) sampling.")
        else:
            lines.append(f"  {prefix}Random window steps (rw_steps = 0): 0 completed random information set steps.")

    if stats and "total_cws" in stats:
        if stats["total_cws"] > 0:
            lines.append(
                f"  {prefix}Codewords accumulated: {stats['total_cws']} distinct "
                f"(w = {stats['min_w']}: {stats['cws']} codewords, {stats['total_hits']} total hits; "
                f"hits min = {stats['hits_min']}, max = {stats['hits_max']}, "
                f"avg = {stats['hits_avg']:.2f}, stdev = {stats['hits_stdev']:.2f})."
            )
            if stats.get("min_hits", 0) > 0 and "avg_hits" in stats:
                import math
                conv_str = "CONVERGED" if stats.get("rw_converged") else "not converged"
                lines.append(
                    f"  {prefix}Hit convergence (min_hits = {stats['min_hits']}, "
                    f"cov_cws = {stats['cov_cws']}): <n> = {stats['avg_hits']:.2f} average hits per codeword "
                    f"({stats['set_cws']} codewords of weight {stats['set_w_lo']}..{stats['set_w_hi']}), "
                    f"P_fail ~ exp(-<n>) = {math.exp(-stats['avg_hits']):.2g} ({conv_str})."
                )
            for ew in stats.get("extra_weights", []):
                lines.append(
                    f"  {prefix}Codewords (w = {ew['w']}): {ew['cws']} codewords, "
                    f"{ew['total_hits']} total hits; hits min = {ew['hits_min']}, "
                    f"max = {ew['hits_max']}, avg = {ew['hits_avg']:.2f}, "
                    f"stdev = {ew['hits_stdev']:.2f}."
                )
            for hw in stats.get("hit_warnings", []):
                lines.append(f"  {prefix}WARNING: {hw}")
        else:
            lines.append(f"  {prefix}Codewords accumulated: 0 distinct non-trivial codewords found.")

    # The RW upper bound is not certified: estimated probability that RW missed a lighter codeword
    if stats and num_rw > 0 and dmax > max(dmin, 1):
        import math
        rw_threads = stats.get("rw_threads", 0)
        step_wall = stats.get("rw_step_time", 0.0) / rw_threads if rw_threads > 0 else 0.0
        avg = stats.get("avg_hits", 0.0)
        if stats.get("min_hits", 0) > 0 and avg > 0.0:
            p_hit = math.exp(-avg)
            if p_hit > 0.01:
                need = math.ceil(num_rw * math.log(100.0) / avg)  # the average number of hits grows with the steps
                t_str = f", ~{_format_duration((need - num_rw) * step_wall)} more" if step_wall > 0.0 else ""
                lines.append(f"  {prefix}Hit-based estimate: P_fail ~ exp(-<n>) = {p_hit:.2g}; 1% after about {need} "
                             f"RW steps in total{t_str}.")
            else:
                lines.append(f"  {prefix}Hit-based estimate: P_fail ~ exp(-<n>) = {p_hit:.2g} (below 1%).")
        n_unif = stats.get("rw_steps_uniform", 0)
        if stats.get("rw_n", 0) > 0 and stats.get("rw_ksub", 0) == 0 and n_unif > 0:
            n_cols, rank = stats["rw_n"], stats["rw_rank"]
            est = rw_infoset_estimate(n_cols, rank, n_unif, dmax, dmin=max(dmin, 1), step_time=step_wall,
                                      steps_total=stats.get("rw_steps", n_unif))
            if est is not None and est["steps_target"] is None:
                lines.append(f"  {prefix}Information-set estimate: a codeword of weight {est['w']} cannot be found "
                             f"by RW (n = {n_cols}, rank(H) = {rank}).")
            elif est is not None:
                t_str = (f", ~{_format_duration(est['time_target'])}" if est["time_target"] is not None else "")
                lines.append(
                    f"  {prefix}Information-set estimate (uniform random information sets; n = {n_cols}, "
                    f"rank(H) = {rank}, {n_unif} of {stats.get('rw_steps', n_unif)} RW steps uniform): a single "
                    f"codeword of weight {est['w']} < {dmax} is found with probability {est['p_find']:.3g} per RW "
                    f"step and would be missed with probability {est['p_miss']:.2g}; 1% after about "
                    f"{est['steps_target']} RW steps in total{t_str}."
                )
                if stats.get("min_w") == dmax and "hits_avg" in stats:
                    n_all = stats.get("rw_steps", n_unif)
                    exp_hits = n_all * rw_infoset_find_prob(n_cols, rank, dmax)
                    verdict = ("the uniform model is conservative here" if stats["hits_avg"] >= exp_hits else
                               "codewords are found less often than in the uniform model, the information-set "
                               "estimate may be optimistic")
                    note = "" if n_unif == n_all else ", including the steps with localized windows"
                    lines.append(
                        f"  {prefix}Model check: expected hits per codeword of weight {dmax} in {n_all} uniform RW "
                        f"steps: {exp_hits:.3g}, observed average {stats['hits_avg']:.2f}{note} ({verdict})."
                    )
        elif stats.get("rw_n", 0) > 0:
            lines.append(f"  {prefix}Information-set estimate: not available (ksub > 0 or kwin > 0: no uniform "
                         f"random information sets).")

    return "\n".join(lines)


def get_cached_distance(
    H: Optional[Any] = None,
    G: Optional[Any] = None,
    L: Optional[Any] = None,
    Hx: Optional[Any] = None,
    Hz: Optional[Any] = None,
    Lx: Optional[Any] = None,
    Lz: Optional[Any] = None,
    dem: Optional[Any] = None,
    circuit: Optional[Any] = None,
    pmin: float = 0.0,
    cache_file: Optional[Union[str, Path]] = None,
    start: Optional[Union[int, str, Sequence[int]]] = None
) -> Optional[Dict[str, Any]]:
    """
    Retrieves the cached distance entry (including bounds and cumulative rw_steps)
    for a given code matrix, CSS code, or DEM.

    If the expert CC start-column list `start` is given (see compute_classical_distance), the
    separate start-list record (cache key suffix ':start=a,b,c') is returned instead of the main record.

    Returns:
        dict with keys {"dist", "dmin", "dmax", "rw_steps", ...} or None if not cached.
    """
    global _distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    if eff_cache_file:
        load_distance_cache(eff_cache_file)
    sfx = _start_key_suffix(_normalize_start_list(start))

    if H is not None:
        if G is not None:
            key = f"quantum:H={get_sparse_array_state(H)}:G={get_sparse_array_state(G)}"
        elif L is not None:
            key = f"quantum:H={get_sparse_array_state(H)}:L={get_sparse_array_state(L)}"
        else:
            key = f"classical:{get_sparse_array_state(H)}"
        entry = _distance_cache.get(key + sfx)
        if entry:
            entry = dict(entry)
            entry["d_info"] = format_bounds_list(entry.get("dmin", 0), entry.get("dmax", 0), entry.get("rw_steps", 0))
        return entry
    elif Hx is not None or Hz is not None:
        hx_st = get_sparse_array_state(Hx) if Hx is not None else "none"
        hz_st = get_sparse_array_state(Hz) if Hz is not None else "none"
        key = f"css:X={hx_st}:Z={hz_st}"
        if Lx is not None or Lz is not None:
            lx_st = get_sparse_array_state(Lx) if Lx is not None else "none"
            lz_st = get_sparse_array_state(Lz) if Lz is not None else "none"
            key = f"{key}:Lx={lx_st}:Lz={lz_st}"
        entry = _distance_cache.get(key + sfx)
        if entry:
            entry = dict(entry)
            if "dmin_X" in entry:
                entry["dX"] = format_bounds_list(
                    entry.get("dmin_X", 0), entry.get("dmax_X", 0), entry.get("rw_steps_X", 0)
                )
            if "dmin_Z" in entry:
                entry["dZ"] = format_bounds_list(
                    entry.get("dmin_Z", 0), entry.get("dmax_Z", 0), entry.get("rw_steps_Z", 0)
                )
        return entry
    elif dem is not None or circuit is not None:
        if dem is None and circuit is not None:
            if hasattr(circuit, 'detector_error_model'):
                obj = circuit.detector_error_model(decompose_errors=True)
            else:
                obj = circuit
        else:
            obj = dem
        dem_st = get_sparse_array_state(obj)
        key = f"dem:{dem_st}" if pmin <= 0.0 else f"dem:{dem_st}:pmin={pmin}"
        entry = _distance_cache.get(key + sfx)
        if entry:
            entry = dict(entry)
            entry["d_info"] = format_bounds_list(entry.get("dmin", 0), entry.get("dmax", 0), entry.get("rw_steps", 0))
        return entry
    return None


# Backward-compatibility aliases
clear_css_distance_cache = clear_distance_cache
enable_css_distance_cache = enable_distance_cache
disable_css_distance_cache = disable_distance_cache


def get_sparse_array_state(A) -> str:
    """Returns a deterministic string representation for JSON-compatible cache keys."""
    if A is None:
        return "none"
    if isinstance(A, (str, Path)):
        path_str = str(Path(A).resolve())
        if os.path.isfile(path_str):
            try:
                with open(path_str, "rb") as f:
                    content_h = hashlib.sha256(f.read()).hexdigest()
                return f"file:{path_str}:{content_h}"
            except Exception:
                return f"file:{path_str}"
        return path_str
    if hasattr(A, 'shape') and hasattr(A, 'dtype') and hasattr(A, 'tobytes'):
        h = hashlib.sha256(A.tobytes()).hexdigest()
        dtype_str = getattr(A.dtype, 'str', str(A.dtype))
        return f"ndarray:{A.shape}:{dtype_str}:{h}"
    if hasattr(A, 'tocsr'):
        csr = A.tocsr()
        h = hashlib.sha256(csr.data.tobytes() + csr.indices.tobytes() + csr.indptr.tobytes()).hexdigest()
        return f"csr:{csr.shape}:{h}"
    if hasattr(A, 'tobytes'):
        h = hashlib.sha256(A.tobytes()).hexdigest()
        return f"bytes:{h}"
    h = hashlib.sha256(str(A).encode('utf-8')).hexdigest()
    return f"str_sha256:{h}"


def create_unique_file(directory: Union[str, Path] = "tmp", extension: str = ".tmp") -> str:
    """Creates a unique temporary file path and ensures the parent directory exists."""
    try:
        os.makedirs(directory, exist_ok=True)
        fd, path = tempfile.mkstemp(suffix=extension, dir=directory)
    except OSError:
        fd, path = tempfile.mkstemp(suffix=extension)
    os.close(fd)
    return path


def _remove_temp_files(files: Sequence[Optional[str]], debug: int = 0) -> None:
    """Removes the temporary files `files` (None entries are skipped); with the debug bit PY_DBG_KEEP_FILES, the
    files are kept and listed instead."""
    for f in files:
        if not f or not os.path.exists(f):
            continue
        if debug & PY_DBG_KEEP_FILES:
            print(f"[dist_m4ri] Kept temporary file: {f}")
            continue
        try:
            os.remove(f)
        except OSError:
            pass


def read_sparse_vectors(filepath: str) -> List[List[int]]:
    """
    Reads a list of sparse vectors from a text file in NZLIST format,
    converting from 1-based indexing (in the file) to 0-based indexing (in Python).

    Args:
        filepath (str): The path to the text file.

    Returns:
        list of list of int: A list where each element is a 0-based sparse vector.
    """
    sparse_vectors = []
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return sparse_vectors

    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            return sparse_vectors
        if first_line != '%% NZLIST':
            raise ValueError(f"Invalid file format in {filepath}: Missing '%% NZLIST' header.")

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            try:
                parts = list(map(int, line.split()))
            except ValueError:
                raise ValueError(f"Non-integer data found on line {line_num}: {line}")

            stated_length = parts[0]
            vector_elements = [x - 1 for x in parts[1:]]
            if len(vector_elements) != stated_length:
                raise ValueError(
                    f"Length mismatch on line {line_num}. "
                    f"Expected {stated_length} elements, but found {len(vector_elements)}."
                )
            sparse_vectors.append(vector_elements)

    return sparse_vectors


def find_dist_m4ri_binary(custom_path: Optional[str] = None) -> str:
    """Finds the dist_m4ri executable."""
    if custom_path:
        if os.path.isfile(custom_path) and os.access(custom_path, os.X_OK):
            return os.path.abspath(custom_path)
        raise FileNotFoundError(
            f"Specified executable '{custom_path}' not found or not executable."
        )

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(pkg_dir, "src", "dist_m4ri"),
        os.path.join(pkg_dir, "dist_m4ri"),
        os.path.join(pkg_dir, "bin", "dist_m4ri"),
        os.path.join(pkg_dir, "..", "dist-m4ri", "src", "dist_m4ri"),
        os.path.join(os.getcwd(), "src", "dist_m4ri"),
        os.path.join(os.getcwd(), "dist_m4ri"),
        os.path.join(os.getcwd(), "bin", "dist_m4ri"),
    ]

    for cand in candidates:
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)

    which_path = shutil.which("dist_m4ri")
    if which_path:
        return which_path

    raise FileNotFoundError(
        "Could not find executable 'dist_m4ri'. Please run 'make -C src' to build it."
    )


def _parse_version(v_str: str) -> Tuple[int, ...]:
    import re
    parts = re.findall(r"\d+", v_str)
    return tuple(int(p) for p in parts) if parts else (0,)


def check_binary_compatibility(binary_path: Optional[str] = None) -> Optional[str]:
    """
    Checks if the backend dist_m4ri binary exists and is compatible (version >= __version__).
    Returns a warning message string if missing or older, or None if compatible (silent).
    """
    try:
        path = find_dist_m4ri_binary(binary_path)
    except (RuntimeError, FileNotFoundError):
        return "Warning: backend binary 'dist_m4ri' not found (run 'make -C src' to build it)"

    try:
        proc = subprocess.run(
            [path, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=2.0
        )
        if proc.returncode == 0 and "version" in proc.stdout:
            bin_ver = proc.stdout.strip().split()[-1]
            if _parse_version(bin_ver) < _parse_version(__version__):
                return (
                    f"Warning: backend binary '{path}' is version {bin_ver} "
                    f"(expected >= {__version__}; run 'make -C src' to rebuild)"
                )
            return None  # Compatible and up-to-date: silent!
        else:
            return (
                f"Warning: backend binary '{path}' does not support --version "
                f"(expected >= {__version__}; run 'make -C src' to rebuild)"
            )
    except Exception as e:
        return f"Warning: failed to check binary '{path}': {e}"


def parse_dist_m4ri_output(stdout: str) -> Tuple[int, int, int]:
    """
    Parses the standard output of dist_m4ri.
    Expected format on stdout: "dmin dmax rw_steps", "dmin dmax", or a single integer.
    
    Returns:
        tuple (dmin, dmax, rw_steps)
    """
    lines = stdout.strip().split('\n')
    for line in reversed(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                return int(parts[0]), int(parts[1]), int(parts[2])
            except ValueError:
                continue
        elif len(parts) == 2:
            try:
                return int(parts[0]), int(parts[1]), 0
            except ValueError:
                continue
        elif len(parts) == 1:
            try:
                val = int(parts[0])
                return val, val, 0
            except ValueError:
                continue

    raise RuntimeError(f"Could not parse dist_m4ri output: {stdout}")


class DistanceResult:
    """
    Structured result for code distance calculations containing:
    - dmin: lower bound on distance
    - dmax: upper bound on distance (minimum non-trivial codeword weight found)
    - rw_steps: cumulative number of completed random window information sets
    - cws: discovered codewords (if requested)
    """
    def __init__(
        self,
        dmin: int,
        dmax: int,
        rw_steps: int = 0,
        cws: Optional[List[List[int]]] = None,
        cws_X: Optional[List[List[int]]] = None,
        cws_Z: Optional[List[List[int]]] = None,
        dmin_X: Optional[int] = None,
        dmax_X: Optional[int] = None,
        rw_steps_X: Optional[int] = None,
        dmin_Z: Optional[int] = None,
        dmax_Z: Optional[int] = None,
        rw_steps_Z: Optional[int] = None,
    ):
        self.dmin = dmin
        self.dmax = dmax
        self.rw_steps = rw_steps
        self.cws = cws
        self.cws_X = cws_X
        self.cws_Z = cws_Z
        self.dmin_X = dmin_X
        self.dmax_X = dmax_X
        self.rw_steps_X = rw_steps_X
        self.dmin_Z = dmin_Z
        self.dmax_Z = dmax_Z
        self.rw_steps_Z = rw_steps_Z

    @property
    def is_exact(self) -> bool:
        return self.dmin > 0 and self.dmin == self.dmax

    @property
    def dist(self) -> int:
        return self.dmin if self.is_exact else (self.dmax if self.dmax > 0 else self.dmin)

    def __int__(self) -> int:
        return self.dist

    def __index__(self) -> int:
        return self.dist

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, DistanceResult):
            return (self.dmin, self.dmax, self.rw_steps) == (other.dmin, other.dmax, other.rw_steps)
        if isinstance(other, (tuple, list)):
            return tuple(self) == tuple(other)
        if isinstance(other, (int, np.integer)):
            return self.dist == other
        return False

    def __iter__(self):
        if self.cws_X is not None or self.cws_Z is not None:
            return iter((self.dmin, self.dmax, self.rw_steps, self.cws_X or [], self.cws_Z or []))
        if self.cws is not None:
            return iter((self.dmin, self.dmax, self.rw_steps, self.cws))
        return iter((self.dmin, self.dmax, self.rw_steps))

    def __getitem__(self, index: int):
        return tuple(self)[index]

    def __len__(self) -> int:
        return len(tuple(self))

    def __str__(self) -> str:
        if self.is_exact:
            return f"{self.dmin} {self.dmax} {self.rw_steps} (exact)"
        return f"{self.dmin} {self.dmax} {self.rw_steps}"

    def __repr__(self) -> str:
        if self.is_exact:
            return f"DistanceResult({self.dmin} {self.dmax} {self.rw_steps} (exact))"
        return f"DistanceResult(dmin={self.dmin}, dmax={self.dmax}, rw_steps={self.rw_steps})"


def check_finc_outc(finC: Optional[str], outC: Optional[str], verbose: bool = False) -> Optional[str]:
    """
    When finC and outC names are identical, an empty or non-existent file is silently ignored
    (with a warning if verbose is True).

    Returns:
        The effective finC filepath to use (or None if ignored).
    """
    if not finC:
        return None
    if outC and (finC == outC or os.path.abspath(finC) == os.path.abspath(outC)):
        if not os.path.exists(finC) or os.path.getsize(finC) == 0:
            if verbose:
                print(f"[dist_m4ri] Warning: finC='{finC}' (identical to outC) is empty or non-existent; "
                      f"silently ignoring input codewords.")
            return None
    return finC


# ---------------------------------------------------------------------------------------------------------------
# Expert CC options: `start` list of CC start columns, `trust_start`, and the disabled options noscan, cbeg, cend.
#
# - start=a,b,c: CC clusters are grown only from the listed columns, and they are NOT limited to larger column
#   indices (unlike cbeg/cend in the dist_m4ri binary).  The lower bound dmin (and an exact result) is valid only
#   if a code symmetry maps every codeword to one whose support contains a listed column (e.g., one column per
#   block of a quasi-cyclic code); the upper bound dmax is always valid.  Results are cached under the separate
#   key '<key>:start=a,b,c' (sorted, distinct columns).  Only dmax and the codewords are merged into the main
#   record, unless the expert option trust_start is set (then dmin and rw_steps are merged as well, and an
#   existing exact start-list record is copied to the main record even if no calculation is needed).
# - noscan, cbeg, cend: disabled in the Python interface (ignored with a warning to stderr); use the dist_m4ri
#   binary directly (noscan: CC only at w=wmax; cbeg/cend: start-column range for split multi-run use, where
#   each column is eventually listed in one of the runs).
# ---------------------------------------------------------------------------------------------------------------

def _normalize_start_list(start: Any) -> Optional[List[int]]:
    """
    Normalizes the expert CC start-column specification `start` into a sorted list of distinct
    non-negative (0-based) column indices, or None if unset.

    Accepts None, an integer, a comma-separated string (e.g., "0,48"), or a sequence of integers.
    A single negative value (e.g., -1, the default of the dist_m4ri binary) means unset (all columns).

    Raises:
        ValueError: on malformed input or on negative entries in a list of several columns.
    """
    import operator
    if start is None:
        return None
    err_msg = f"invalid start={start!r}: expected a comma-separated list of column indices, e.g., start=0,48"
    vals: List[int] = []
    if isinstance(start, str):
        s = start.strip()
        if not s:
            return None
        try:
            vals = [int(tok) for tok in s.split(",")]
        except ValueError:
            raise ValueError(err_msg) from None
    else:
        try:
            items = list(start)
        except TypeError:
            items = [start]
        for v in items:
            if isinstance(v, (bool, str, bytes)):
                raise ValueError(err_msg)
            try:
                vals.append(operator.index(v))
            except TypeError:
                raise ValueError(err_msg) from None
    if not vals:
        return None
    if len(vals) == 1 and vals[0] < 0:
        return None
    if any(v < 0 for v in vals):
        raise ValueError(f"invalid start={start!r}: column indices must be non-negative")
    return sorted(set(vals))


def _start_key_suffix(start_list: Optional[List[int]]) -> str:
    """Returns the cache key suffix ':start=a,b,c' for a normalized start list ('' if unset)."""
    if not start_list:
        return ""
    return ":start=" + ",".join(str(c) for c in start_list)


def _format_start_list(start_list: List[int], max_show: int = 8) -> str:
    """Formats a start list for messages, showing at most `max_show` entries (as the dist_m4ri binary)."""
    shown = ",".join(str(c) for c in start_list[:max_show])
    return shown + (",..." if len(start_list) > max_show else "")


def _warn_disabled_options(noscan: Any = 0, cbeg: Any = None, cend: Any = None) -> None:
    """
    Issues a warning to stderr for each of the expert options noscan, cbeg, and cend that is set.
    These options are disabled in the Python interface: they are accepted for backward compatibility
    but ignored, since their results are not valid distance bounds on their own.
    """
    for name, val, is_flag in (("noscan", noscan, True), ("cbeg", cbeg, False), ("cend", cend, False)):
        if val is None:
            continue
        try:
            ival = int(val)
            active = (ival != 0) if is_flag else (ival >= 0)
        except (TypeError, ValueError):
            active = True
        if active:
            sys.stderr.write(
                f"# Warning: {name}={val} is disabled in the Python interface (ignored); "
                f"use the dist_m4ri binary directly\n"
            )


def _start_warning_text(start_list: List[int]) -> str:
    """Returns the stderr warning for a CC start list (same wording as the dist_m4ri binary)."""
    return (
        f"# WARNING: start={_format_start_list(start_list)} (expert option): CC clusters are grown only from "
        f"{len(start_list)} listed column(s)\n"
        f"#   (not limited to larger column indices); dmin is valid only if a code symmetry maps every\n"
        f"#   minimum-weight codeword to one containing a listed column\n"
    )


def _ksub_warning_text(ksub: int) -> str:
    """Returns the stderr warning for the experimental RW option ksub (same wording as the dist_m4ri binary)."""
    return (
        f"# WARNING: ksub={ksub} is experimental and should not be used: a RW step reduces the span of ksub rows\n"
        f"#   of a fixed basis of ker(H), so that a few codewords are found much more often than others; with\n"
        f"#   min_hits, RW may end early, and exp(-<n>) underestimates the probability to miss a lighter codeword\n"
    )


def _warn_experimental_ksub(ksub: Any, method: int, solver: str = "dist_m4ri") -> None:
    """
    Issues the warning for the experimental RW option ksub > 0 to stderr (with or without min_hits), if RW runs
    (method=1 or 3) with the dist_m4ri binary.
    """
    try:
        ival = int(ksub or 0)
    except (TypeError, ValueError):
        return
    if ival > 0 and (int(method) & 1) and solver != "codedistance":
        sys.stderr.write(_ksub_warning_text(ival))


def _prepare_start_option(
    start: Any,
    trust_start: bool = False,
    noscan: Any = 0,
    cbeg: Any = None,
    cend: Any = None,
    solver: str = "dist_m4ri"
) -> Optional[List[int]]:
    """
    Validates the expert CC options of a compute_*_distance() call and issues the warnings to stderr.

    Returns:
        The normalized start list, or None if unset (or not applicable to the selected solver).
    """
    _warn_disabled_options(noscan, cbeg, cend)
    start_list = _normalize_start_list(start)
    if start_list is None:
        if trust_start:
            sys.stderr.write("# Warning: trust_start=1 has no effect without a start list\n")
        return None
    if solver == "codedistance":
        sys.stderr.write(f"# Warning: start={_format_start_list(start_list)} is ignored with solver='codedistance'\n")
        return None
    msg = _start_warning_text(start_list)
    if _use_distance_cache and trust_start:
        msg += "#   trust_start=1: start-list results are accepted as valid and copied to the main cache record\n"
    elif _use_distance_cache:
        msg += (
            "#   results are cached under a separate ':start=...' key; only dmax and codewords update the\n"
            "#   main cache record (use trust_start=1 to also accept dmin)\n"
        )
    sys.stderr.write(msg)
    return start_list


def _seed_bounds_from_entry(
    eff_dmin: int, eff_dmax: int, entry: Optional[Dict[str, Any]], sfx: str = ""
) -> Tuple[int, int]:
    """
    Tightens the bounds (eff_dmin, eff_dmax) with those stored in a cache entry (for a CSS entry, use the
    sector suffix sfx='_X' or '_Z'; missing sector fields fall back to the combined 'dmin' / 'dmax').
    """
    if not entry:
        return eff_dmin, eff_dmax
    c_dmin = entry.get(f"dmin{sfx}", entry.get("dmin", 0)) or 0
    c_dmax = entry.get(f"dmax{sfx}", entry.get("dmax", 0)) or 0
    if c_dmax > 0:
        eff_dmax = min(eff_dmax, c_dmax) if eff_dmax > 0 else c_dmax
    if c_dmin > 1:
        eff_dmin = max(eff_dmin, c_dmin) if eff_dmin > 1 else c_dmin
    return eff_dmin, eff_dmax


def _union_cws(base: Optional[List[List[int]]], extra: Optional[List[List[int]]]) -> List[List[int]]:
    """Returns the union of two codeword lists (duplicates removed, sorted by weight)."""
    combined = [list(cw) for cw in (base or [])]
    seen = {tuple(cw) for cw in combined}
    for cw in (extra or []):
        t = tuple(cw)
        if t not in seen:
            combined.append(list(cw))
            seen.add(t)
    combined.sort(key=len)
    return combined


def _single_entry_is_exact(entry: Dict[str, Any], need_cws: bool = False) -> bool:
    """Cache hit test for a classical / quantum / DEM entry (exact distance, and codewords if needed)."""
    if not (entry.get("dmin", 0) > 0 and entry.get("dmin") == entry.get("dmax")):
        return False
    return (not need_cws) or bool(entry.get("cws"))


def _css_entry_is_exact(entry: Dict[str, Any], can_x: bool, can_z: bool, need_cws: bool = False) -> bool:
    """Cache hit test for a CSS entry (exact distance in all computed sectors, and codewords if needed)."""
    c_dmin_x = entry.get("dmin_X", entry.get("dmin", 0))
    c_dmax_x = entry.get("dmax_X", entry.get("dmax", 0))
    c_dmin_z = entry.get("dmin_Z", entry.get("dmin", 0))
    c_dmax_z = entry.get("dmax_Z", entry.get("dmax", 0))
    x_exact = (not can_x) or (c_dmin_x > 0 and c_dmin_x == c_dmax_x)
    z_exact = (not can_z) or (c_dmin_z > 0 and c_dmin_z == c_dmax_z)
    if not (x_exact and z_exact and (c_dmax_x > 0 or c_dmax_z > 0)):
        return False
    return (not need_cws) or bool(entry.get("cws_X") and entry.get("cws_Z"))


def _merge_start_into_plain(
    plain_key: str,
    src: Dict[str, Any],
    trust: bool = False,
    rw_delta: Optional[Dict[str, int]] = None,
    sectors: Optional[Sequence[str]] = None
) -> bool:
    """
    Merges the start-list cache record `src` (key plain_key + ':start=...') into the main record `plain_key`.

    The upper bound dmax (minimum of the positive values) and the codewords (union) are always valid and
    are always merged.  Only with `trust` (expert option trust_start) the lower bound dmin is merged as well
    (maximum, capped at dmax), and rw_steps are either incremented by the steps of the run, `rw_delta`
    (after a calculation), or set to the maximum of the two records (copy of a cached record, rw_delta=None).

    Args:
        plain_key: Main cache key.
        src: Start-list cache record.
        trust: Accept the start-list lower bound as valid.
        rw_delta: RW steps of the run by key suffix ('' for a single record; '_X', '_Z' for CSS sectors).
        sectors: None for a single (classical, quantum, or DEM) record; for a CSS record, the list of
            computed sectors among 'X' and 'Z'.

    Returns:
        True if the main record was created or modified.
    """
    global _distance_cache
    prev = _distance_cache.get(plain_key)
    new: Dict[str, Any] = dict(prev) if prev else {}
    sfx_list = [""] if sectors is None else [f"_{s}" for s in sectors]
    useful = trust
    res: Dict[str, Tuple[int, int, int]] = {}
    for sfx in sfx_list:
        s_dmin = src.get(f"dmin{sfx}", 0) or 0
        s_dmax = src.get(f"dmax{sfx}", 0) or 0
        p_dmin = new.get(f"dmin{sfx}", 0) or 0
        p_dmax = new.get(f"dmax{sfx}", 0) or 0
        p_rw = new.get(f"rw_steps{sfx}", 0) or 0
        dmax = min(p_dmax, s_dmax) if (p_dmax > 0 and s_dmax > 0) else (p_dmax if p_dmax > 0 else s_dmax)
        dmin, rw = p_dmin, p_rw
        if trust:
            dmin = max(p_dmin, s_dmin)
            if rw_delta is not None:
                rw = p_rw + (rw_delta.get(sfx, 0) or 0)
            else:
                rw = max(p_rw, src.get(f"rw_steps{sfx}", 0) or 0)
        if dmax > 0 and dmin > dmax:
            dmin = dmax
        if sectors is not None and dmin > 0 and dmin == dmax:
            rw = 0  # CSS sector convention (see compute_css_distance)
        s_cws = src.get(f"cws{sfx}", []) or []
        if s_dmax > 0 or s_cws:
            useful = True
        new[f"cws{sfx}"] = _union_cws(new.get(f"cws{sfx}", []), s_cws)
        res[sfx] = (dmin, dmax, rw)
    if not prev and not useful:
        return False

    if sectors is None:
        dmin, dmax, rw = res[""]
        new["dist"] = dmin if (dmin == dmax or dmax == 0) else dmax
        new["dmin"], new["dmax"], new["rw_steps"] = dmin, dmax, rw
        new["d_info"] = format_bounds_list(dmin, dmax, rw)
    else:
        for s in ("X", "Z"):
            sfx = f"_{s}"
            if sfx in res:
                dmin, dmax, rw = res[sfx]
                new[f"dmin{sfx}"], new[f"dmax{sfx}"], new[f"rw_steps{sfx}"] = dmin, dmax, rw
                new[f"d{s}"] = format_bounds_list(dmin, dmax, rw)
            else:
                for fld, dflt in (("dmin", 0), ("dmax", 0), ("rw_steps", 0), ("cws", [])):
                    new.setdefault(f"{fld}{sfx}", dflt)
                new.setdefault(f"d{s}", None)
        vals = list(res.values())
        dists = [dmin if (dmin == dmax or dmax == 0) else dmax for dmin, dmax, _ in vals]
        new["dist"] = min(dists) if all(d > 0 for d in dists) else max(dists)
        new["dmin"] = min(v[0] for v in vals)
        new["dmax"] = min(v[1] for v in vals) if all(v[1] > 0 for v in vals) else 0
        new["rw_steps"] = sum(v[2] for v in vals)
    _distance_cache[plain_key] = new
    return True


def _resolve_start_cache(
    plain_key: str,
    start_list: Optional[List[int]],
    is_hit: Callable[[Dict[str, Any]], bool],
    trust_start: bool,
    merge_copy: Callable[[Dict[str, Any]], bool],
    eff_cache_file: Optional[str],
    verbose: bool = False
) -> Tuple[str, Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """
    Resolves the cache record of a run with an optional expert CC start list.

    Without a start list, or if the main record already holds the requested exact result, the main record
    is used.  Otherwise the separate start-list record (key plain_key + ':start=a,b,c') is used; with
    trust_start, an existing exact start-list record is copied to the main record (via `merge_copy`), and
    the cache file is saved, even though no calculation is needed.

    Returns:
        (code_key, cached_entry, plain_entry): the key and the existing record (or None) for the results of
        this run, and the main record (or None) used to seed the bounds of a start-list run.
    """
    plain_entry = _distance_cache.get(plain_key)
    if not start_list:
        return plain_key, plain_entry, None
    if plain_entry is not None and is_hit(plain_entry):
        if verbose:
            print(f"[dist_m4ri] Start list not needed: the main record '{plain_key}' has the exact distance")
        return plain_key, plain_entry, None
    code_key = plain_key + _start_key_suffix(start_list)
    cached_entry = _distance_cache.get(code_key)
    if trust_start and cached_entry is not None and is_hit(cached_entry):
        if merge_copy(cached_entry) and eff_cache_file:
            save_distance_cache(eff_cache_file)
        if verbose:
            print(f"[dist_m4ri] trust_start: copied the exact start-list record to the main record '{plain_key}'")
    return code_key, cached_entry, plain_entry


def run_dist_m4ri(
    dist_m4ri_path: Optional[str] = None,
    method: int = 3,
    finH: Optional[str] = None,
    finG: Optional[str] = None,
    finL: Optional[str] = None,
    fin: Optional[str] = None,
    finC: Optional[str] = None,
    fdem: Optional[str] = None,
    dmin: int = 0,
    dmax: int = 0,
    wmax: int = 0,
    wmin: int = 1,
    dexp: int = 0,
    dest: int = 0,
    steps: Optional[int] = None,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    smax: Optional[int] = None,
    start: Optional[Union[int, str, Sequence[int]]] = None,
    cbeg: Optional[int] = None,
    cend: Optional[int] = None,
    css: Optional[int] = None,
    noscan: int = 0,
    classical: int = -1,
    dW: int = -1,
    maxC: int = 0,
    pmin: float = 0.0,
    outC: Optional[str] = None,
    seed: int = 0,
    debug: int = 0,
    nothrottle: bool = False,
    chunk_size: int = 0,
    ksub: int = 0,
    kwin: int = 0,
    win_mode: int = 0,
    min_hits: Optional[int] = None,
    cov_cws: int = 100,
    refresh: int = 0,
    verbose: bool = False,
    stop_event: Optional[threading.Event] = None,
    warn_start: bool = True,
    warn_ksub: bool = True
) -> Tuple[int, int, int]:
    """
    Low-level invocation of the multithreaded dist_m4ri binary.

    Expert CC options:
        start: List of CC start columns (int, comma-separated string "a,b,c", or sequence of ints;
            None or a single negative value: all columns), passed as start=a,b,c (sorted, distinct).
            CC clusters are grown only from the listed columns, and they are not limited to larger
            column indices.  WARNING: the returned dmin (and an exact result) is valid only if a code
            symmetry maps every minimum-weight codeword to one containing a listed column (e.g., one
            column per block of a quasi-cyclic code); dmax is always valid.
        warn_start: If True (default), issue the start-list warning to stderr.
        noscan / cbeg / cend: Disabled in the Python interface: accepted for backward compatibility,
            but ignored with a warning to stderr (use the dist_m4ri binary directly).

    Experimental RW option:
        ksub: Experimental, should not be used: subspace dimension sampled from ker(H) in RW (0: full-matrix RW).
            A RW step reduces the span of ksub rows of a fixed basis of ker(H), so that a few codewords are found
            much more often than others (min_hits may end RW early).
        warn_ksub: If True (default), issue the warning for ksub > 0 to stderr (with method=1 or 3).

    Debug output:
        debug: Debug bitmap; the bits in BIN_DBG_MASK are passed to the binary (its diagnostic output on stderr
            is printed), and PY_DBG_COMMANDS echoes the command line and the run time.  verbose=True implies
            the bits VERBOSE_DEBUG.

    Returns:
        tuple (dmin, dmax, rw_steps)
    """
    global _last_run_stats
    exec_path = find_dist_m4ri_binary(dist_m4ri_path)

    start_list = _normalize_start_list(start)
    _warn_disabled_options(noscan, cbeg, cend)
    if start_list and warn_start:
        sys.stderr.write(_start_warning_text(start_list))
    if warn_ksub:
        _warn_experimental_ksub(ksub, method)

    finC = check_finc_outc(finC, outC, verbose=False)

    if method == 2 and wmax <= 0:
        if dmax > 0:
            wmax = dmax
        elif timeout <= 0.0:
            raise ValueError("either parameter wmax>0, dmax>0, or timeout>0 should be specified for CC method=2.")

    eff_debug = (debug | VERBOSE_DEBUG) if verbose else debug
    bin_debug = eff_debug & BIN_DBG_MASK  # the higher bits are used by the Python wrapper only
    cmd = [exec_path, f"debug={bin_debug}", f"method={method}"]

    if finH: cmd.append(f"finH={finH}")
    if finG: cmd.append(f"finG={finG}")
    if finL: cmd.append(f"finL={finL}")
    if fin: cmd.append(f"fin={fin}")
    if finC: cmd.append(f"finC={finC}")
    if fdem: cmd.append(f"fdem={fdem}")
    if dmin > 0: cmd.append(f"dmin={dmin}")
    if dmax > 0: cmd.append(f"dmax={dmax}")
    if wmax > 0: cmd.append(f"wmax={wmax}")
    if wmin is not None and wmin != 1: cmd.append(f"wmin={wmin}")
    if dexp > 0: cmd.append(f"dexp={dexp}")
    elif dest > 0: cmd.append(f"dest={dest}")
    if steps is not None and steps >= 0: cmd.append(f"steps={steps}")
    if threads is not None and threads > 0: cmd.append(f"threads={threads}")
    if timeout is not None and timeout >= 0: cmd.append(f"timeout={timeout}")
    if smax is not None: cmd.append(f"smax={smax}")
    if start_list: cmd.append("start=" + ",".join(str(c) for c in start_list))
    if css is not None: cmd.append(f"css={css}")
    if classical >= 0: cmd.append(f"classical={classical}")
    if dW >= 0: cmd.append(f"dW={dW}")
    if maxC > 0: cmd.append(f"maxC={maxC}")
    if pmin > 0.0: cmd.append(f"pmin={pmin}")
    if outC: cmd.append(f"outC={outC}")
    if seed != 0: cmd.append(f"seed={seed}")
    if nothrottle: cmd.append("nothrottle=1")
    if chunk_size > 0: cmd.append(f"chunk_size={chunk_size}")
    if ksub > 0: cmd.append(f"ksub={ksub}")
    if kwin > 0: cmd.append(f"kwin={kwin}")
    if win_mode != 0: cmd.append(f"win_mode={win_mode}")
    if min_hits is not None and min_hits >= 0: cmd.append(f"min_hits={min_hits}")
    if cov_cws != 100: cmd.append(f"cov_cws={cov_cws}")
    if refresh > 0: cmd.append(f"refresh={refresh}")

    if eff_debug & PY_DBG_COMMANDS:
        print(f"[dist_m4ri] Running: {' '.join(cmd)}")

    t_start = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if stop_event is not None:
        while proc.poll() is None:
            if stop_event.is_set():
                proc.terminate()
                try:
                    proc.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    proc.kill()
                raise RuntimeError("dist_m4ri execution cancelled by stop_event")
            time.sleep(0.05)
        stdout, stderr = proc.communicate()
    else:
        stdout, stderr = proc.communicate()

    if eff_debug & PY_DBG_COMMANDS:
        print(f"[dist_m4ri] Finished in {time.time() - t_start:.3f} s (exit code {proc.returncode})")

    if proc.returncode != 0:
        raise RuntimeError(f"dist_m4ri failed with exit code {proc.returncode}:\n{stderr}")

    _last_run_stats = _parse_stderr_stats(stderr)
    if (verbose or bin_debug > 0) and stderr:
        print(stderr.rstrip())

    return parse_dist_m4ri_output(stdout)


def _matrix_to_file(matrix, extension: str = ".mtx", temp_dir: str = "tmp") -> str:
    """Helper to convert a matrix (numpy or scipy sparse) or file path to an MTX file path."""
    if isinstance(matrix, (str, Path)):
        return str(matrix)

    import numpy as np
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix, issparse

    path = create_unique_file(directory=temp_dir, extension=extension)
    if issparse(matrix):
        csr = matrix.astype(np.int8)
        mmwrite(path, csr, symmetry='general')
    else:
        mat_arr = np.asarray(matrix, dtype=np.int8)
        csr = csr_matrix(mat_arr)
        mmwrite(path, csr, symmetry='general')
    return path


def has_noise(circuit: Any) -> bool:
    """Checks if a Stim circuit contains any noise instructions."""
    try:
        stim = _get_stim()
    except ImportError:
        return False
    if not isinstance(circuit, stim.Circuit):
        return False
    noisy_gate_names = {
        "DEPOLARIZE1", "DEPOLARIZE2", "PAULI_CHANNEL_1", "PAULI_CHANNEL_2",
        "X_ERROR", "Y_ERROR", "Z_ERROR", "E", "ELSE_CORRELATED_ERROR"
    }
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            if has_noise(inst.body_copy()):
                return True
        elif inst.name in noisy_gate_names:
            return True
    return False


def _add_noise_recursive(
    circuit: Any, p: float, num_qubits: int, active_qubits: Set[int]
) -> Any:
    """Recursively adds uniform circuit-level noise to a Stim circuit (see add_noise())."""
    stim = _get_stim()
    noisy_circuit = stim.Circuit()
    annotations = {
        "QUBIT_COORDS", "DETECTOR", "OBSERVABLE_INCLUDE", "TICK", "SHIFT_COORDS"
    }
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            noisy_body = _add_noise_recursive(
                inst.body_copy(), p, num_qubits, active_qubits
            )
            noisy_circuit.append(
                stim.CircuitRepeatBlock(inst.repeat_count, noisy_body)
            )
        elif inst.name == "TICK":
            idle_qubits = set(range(num_qubits)) - active_qubits
            if idle_qubits:
                noisy_circuit.append(
                    "DEPOLARIZE1", sorted(list(idle_qubits)), 0.1 * p
                )
            noisy_circuit.append(inst)
            active_qubits.clear()
        else:
            if inst.name not in annotations:
                for t in inst.targets_copy():
                    if t.value >= 0:
                        active_qubits.add(t.value)

            qubit_targets = [
                t.value for t in inst.targets_copy() if t.value >= 0
            ]
            if inst.name == "RX":
                noisy_circuit.append(inst)
                if qubit_targets:
                    noisy_circuit.append("Z_ERROR", qubit_targets, p)
            elif inst.name in ["R", "RZ", "RY"]:
                noisy_circuit.append(inst)
                if qubit_targets:
                    noisy_circuit.append("X_ERROR", qubit_targets, p)
            elif inst.name == "MX":
                if qubit_targets:
                    noisy_circuit.append("Z_ERROR", qubit_targets, p)
                noisy_circuit.append(inst)
            elif inst.name in ["M", "MZ", "MY"]:
                if qubit_targets:
                    noisy_circuit.append("X_ERROR", qubit_targets, p)
                noisy_circuit.append(inst)
            elif inst.name == "MRX":
                if qubit_targets:
                    noisy_circuit.append("Z_ERROR", qubit_targets, p)
                noisy_circuit.append(inst)
                if qubit_targets:
                    noisy_circuit.append("Z_ERROR", qubit_targets, p)
            elif inst.name in ["MR", "MRZ", "MRY"]:
                if qubit_targets:
                    noisy_circuit.append("X_ERROR", qubit_targets, p)
                noisy_circuit.append(inst)
                if qubit_targets:
                    noisy_circuit.append("X_ERROR", qubit_targets, p)
            else:
                noisy_circuit.append(inst)
                if inst.name in [
                    "I", "X", "Y", "Z", "H", "S", "S_DAG",
                    "SQRT_X", "SQRT_X_DAG", "SQRT_Y", "SQRT_Y_DAG",
                    "SQRT_Z", "SQRT_Z_DAG"
                ]:
                    noisy_circuit.append(
                        "DEPOLARIZE1", inst.targets_copy(), 0.1 * p
                    )
                elif inst.name in [
                    "CX", "CY", "CZ", "SWAP", "XCZ", "YCX", "YCY", "YCZ"
                ]:
                    noisy_circuit.append("DEPOLARIZE2", inst.targets_copy(), p)
    return noisy_circuit


def add_noise(circuit: Any, p: float = 0.001) -> Any:
    """
    Adds uniform circuit-level noise to a Stim circuit: DEPOLARIZE2(p) after two-qubit gates, DEPOLARIZE1(p/10) after
    single-qubit gates and on the idle qubits in each TICK, and X (Z) flips with probability p after Z-basis (X-basis)
    resets and before measurements.

    Args:
        circuit: stim.Circuit or path to a .stim file.
        p: Noise strength (for the distance, only which errors are possible matters).

    Returns:
        The noisy stim.Circuit.
    """
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))
    num_qubits = circuit.num_qubits
    active_qubits: Set[int] = set()
    noisy_circuit = _add_noise_recursive(circuit, p, num_qubits, active_qubits)
    if active_qubits:
        final_idle_qubits = set(range(num_qubits)) - active_qubits
        if final_idle_qubits:
            noisy_circuit.append(
                "DEPOLARIZE1", sorted(list(final_idle_qubits)), p
            )
    return noisy_circuit


def remove_empty_detectors(circuit: Any) -> Tuple[Any, int]:
    """Removes DETECTOR instructions that have no measurement record targets."""
    stim = _get_stim()
    new_circuit = stim.Circuit()
    removed = 0
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            body, r = remove_empty_detectors(inst.body_copy())
            removed += r * inst.repeat_count
            new_circuit.append(stim.CircuitRepeatBlock(inst.repeat_count, body))
        elif inst.name == "DETECTOR":
            rec_targets = [
                t for t in inst.targets_copy()
                if t.is_measurement_record_target
            ]
            if not rec_targets:
                removed += 1
            else:
                new_circuit.append(inst)
        else:
            new_circuit.append(inst)
    return new_circuit, removed


def has_repeat_block(circuit: Any) -> bool:
    """Checks if a Stim circuit contains at least one REPEAT block."""
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))
    if not isinstance(circuit, stim.Circuit):
        return False
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            return True
    return False


def set_circuit_rounds(circuit: Any, rounds: int) -> Any:
    """Sets the repeat_count of CircuitRepeatBlocks in a Stim circuit.

    For circuits with a single repeat block, sets its repeat_count to rounds.
    For circuits with two repeat blocks (e.g. preamble + periodic body),
    adjusts them to achieve total requested rounds.

    Args:
        circuit: A stim.Circuit instance.
        rounds: The new repetition count to set for repeat blocks. If <= 0,
            repeat blocks are omitted.

    Returns:
        A new stim.Circuit instance with modified repeat counts.
    """
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))
    if not isinstance(circuit, stim.Circuit):
        return circuit
    if rounds < 0:
        raise ValueError(f"Invalid rounds={rounds}; must be >= 0.")

    top_blocks = [
        inst.repeat_count for inst in circuit
        if isinstance(inst, stim.CircuitRepeatBlock)
    ]
    if len(top_blocks) > 1:
        sys.stderr.write(
            f"# Warning: Circuit has {len(top_blocks)} repeat blocks "
            f"({top_blocks}) which may not be handled correctly by this "
            "script.\n"
        )

    new_circuit = stim.Circuit()
    if len(top_blocks) == 2:
        preamble = top_blocks[0]
        b1_count = min(preamble, rounds)
        b2_count = max(0, rounds - preamble)
        block_idx = 0
        for inst in circuit:
            if isinstance(inst, stim.CircuitRepeatBlock):
                cnt = b1_count if block_idx == 0 else b2_count
                block_idx += 1
                if cnt > 0:
                    new_body = set_circuit_rounds(inst.body_copy(), cnt)
                    new_circuit.append(stim.CircuitRepeatBlock(cnt, new_body))
            else:
                new_circuit.append(inst)
    else:
        for inst in circuit:
            if isinstance(inst, stim.CircuitRepeatBlock):
                if rounds > 0:
                    new_body = set_circuit_rounds(inst.body_copy(), rounds)
                    new_circuit.append(
                        stim.CircuitRepeatBlock(rounds, new_body)
                    )
            else:
                new_circuit.append(inst)
    return new_circuit


def count_circuit_rounds(circuit: Any) -> int:
    """Returns total repeat count of all repeat blocks in circuit, or 1 if unrolled."""
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))
    if not isinstance(circuit, stim.Circuit):
        return 1
    blocks = [
        inst.repeat_count for inst in circuit
        if isinstance(inst, stim.CircuitRepeatBlock)
    ]
    if len(blocks) > 1:
        sys.stderr.write(
            f"# Warning: Circuit has {len(blocks)} repeat blocks "
            f"({blocks}) which may not be handled correctly by this "
            "script.\n"
        )
    return sum(blocks) if blocks else 1


def detect_basis(
    circuit: Any,
    filepath: Optional[Union[str, Path]] = None,
    basis_arg: Optional[str] = None
) -> str:
    """Determines the primary memory basis ('X' or 'Z') of a Stim circuit."""
    if basis_arg is not None and str(basis_arg).upper() in ["X", "Z"]:
        return str(basis_arg).upper()

    if filepath is not None:
        stem = Path(filepath).stem.upper()
        parts = stem.replace("-", "_").split("_")
        if "X" in parts or "HX" in parts or stem.endswith("X"):
            return "X"
        if "Z" in parts or "HZ" in parts or stem.endswith("Z"):
            return "Z"

    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))

    rx_count = 0
    rz_count = 0
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            break
        if inst.name == "RX":
            rx_count += len(inst.targets_copy())
        elif inst.name in ["R", "RZ"]:
            rz_count += len(inst.targets_copy())
        elif inst.name in ["M", "MX", "MZ", "MR", "MRX", "MRZ"]:
            break
    if rx_count > rz_count:
        return "X"
    if rz_count > rx_count:
        return "Z"

    mx_count = 0
    mz_count = 0
    for inst in reversed(list(circuit)):
        if isinstance(inst, stim.CircuitRepeatBlock):
            break
        if inst.name in ["MX", "MRX"]:
            mx_count += len(inst.targets_copy())
        elif inst.name in ["M", "MZ", "MR", "MRZ"]:
            mz_count += len(inst.targets_copy())
    if mx_count > mz_count:
        return "X"
    return "Z"


def track_circuit_carriers(circuit: Any) -> List[List[int]]:
    """Dynamically tracks logical syndrome carriers across the circuit to measurements.

    Accounts for physical qubit permutations via SWAP, CXSWAP, CZSWAP, ISWAP,
    and ISWAP_DAG gates, as well as CNOT/CZ syndrome handoffs onto fresh ancilla
    or routing qubits.

    Returns:
        meas_carriers: List mapping each measurement index to a list of initial
            logical qubit IDs whose syndrome or state is carried into that
            measurement.
    """
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))
    num_qubits = circuit.num_qubits
    phys_to_log = {q: q for q in range(num_qubits)}

    tot_resets: Dict[int, int] = {q: 0 for q in range(num_qubits)}
    tot_meas: Dict[int, int] = {q: 0 for q in range(num_qubits)}

    def pass1(block: Any, p_to_l: Dict[int, int], mult: int = 1) -> None:
        for inst in block:
            if isinstance(inst, stim.CircuitRepeatBlock):
                rep_count = inst.repeat_count
                body = inst.body_copy()
                map_copy = dict(p_to_l)
                for sub_inst in body:
                    if not isinstance(sub_inst, stim.CircuitRepeatBlock):
                        if sub_inst.name in [
                            "SWAP", "CXSWAP", "CZSWAP", "ISWAP", "ISWAP_DAG"
                        ]:
                            t_vals = [
                                t.value for t in sub_inst.targets_copy()
                                if t.is_qubit_target
                            ]
                            for idx in range(0, len(t_vals), 2):
                                q1, q2 = t_vals[idx], t_vals[idx + 1]
                                map_copy[q1], map_copy[q2] = (
                                    map_copy[q2], map_copy[q1]
                                )
                if map_copy == p_to_l:
                    pass1(body, p_to_l, mult * rep_count)
                else:
                    for _ in range(rep_count):
                        pass1(body, p_to_l, mult)
                continue
            name = inst.name
            t_vals = [
                t.value for t in inst.targets_copy() if t.is_qubit_target
            ]
            if name in [
                "R", "RX", "RY", "RZ", "MR", "MRX", "MRY", "MRZ"
            ]:
                for q in t_vals:
                    tot_resets[p_to_l[q]] += mult
            if name in [
                "M", "MX", "MY", "MZ", "MR", "MRX", "MRY", "MRZ"
            ]:
                for q in t_vals:
                    tot_meas[p_to_l[q]] += mult
            if name in ["SWAP", "CXSWAP", "CZSWAP", "ISWAP", "ISWAP_DAG"]:
                for idx in range(0, len(t_vals), 2):
                    q1, q2 = t_vals[idx], t_vals[idx + 1]
                    p_to_l[q1], p_to_l[q2] = p_to_l[q2], p_to_l[q1]

    pass1(circuit, dict(phys_to_log), 1)
    max_resets = max(tot_resets.values()) if tot_resets else 0
    max_meas = max(tot_meas.values()) if tot_meas else 0
    data_log = set()
    for q in range(num_qubits):
        if max_resets > 1:
            if tot_resets[q] <= 1:
                data_log.add(q)
        elif max_meas > 1:
            if tot_meas[q] <= 1:
                data_log.add(q)

    active_carriers: Dict[int, Set[int]] = {
        q: {q} for q in range(num_qubits)
    }
    fresh_anc: Set[int] = set(range(num_qubits)) - data_log
    meas_carriers: List[List[int]] = []

    def pass2(block: Any, p_to_l: Dict[int, int]) -> None:
        for inst in block:
            if isinstance(inst, stim.CircuitRepeatBlock):
                rep_count = inst.repeat_count
                body = inst.body_copy()
                if rep_count <= 0:
                    continue
                map_copy = dict(p_to_l)
                for sub_inst in body:
                    if not isinstance(sub_inst, stim.CircuitRepeatBlock):
                        if sub_inst.name in [
                            "SWAP", "CXSWAP", "CZSWAP", "ISWAP", "ISWAP_DAG"
                        ]:
                            t_vals = [
                                t.value for t in sub_inst.targets_copy()
                                if t.is_qubit_target
                            ]
                            for idx in range(0, len(t_vals), 2):
                                q1, q2 = t_vals[idx], t_vals[idx + 1]
                                map_copy[q1], map_copy[q2] = (
                                    map_copy[q2], map_copy[q1]
                                )
                if map_copy == p_to_l:
                    idx0 = len(meas_carriers)
                    pass2(body, p_to_l)
                    idx1 = len(meas_carriers)
                    if rep_count > 1 and idx1 > idx0:
                        body_carriers = meas_carriers[idx0:idx1]
                        meas_carriers.extend(body_carriers * (rep_count - 1))
                else:
                    for _ in range(rep_count):
                        pass2(body, p_to_l)
                continue
            name = inst.name
            t_vals = [
                t.value for t in inst.targets_copy() if t.is_qubit_target
            ]
            if name in ["R", "RX", "RY", "RZ"]:
                for q in t_vals:
                    lq = p_to_l[q]
                    active_carriers[lq] = {lq}
                    if lq not in data_log:
                        fresh_anc.add(lq)
            elif name in ["M", "MX", "MY", "MZ"]:
                for q in t_vals:
                    lq = p_to_l[q]
                    meas_carriers.append(sorted(list(active_carriers[lq])))
            elif name in ["MR", "MRX", "MRY", "MRZ"]:
                for q in t_vals:
                    lq = p_to_l[q]
                    meas_carriers.append(sorted(list(active_carriers[lq])))
                    active_carriers[lq] = {lq}
                    if lq not in data_log:
                        fresh_anc.add(lq)
            elif name in [
                "CX", "ZCX", "CY", "ZCY", "CZ", "ZCZ",
                "XCX", "XCY", "XCZ", "YCX", "YCY", "YCZ"
            ]:
                for idx in range(0, len(t_vals), 2):
                    q1, q2 = t_vals[idx], t_vals[idx + 1]
                    l1, l2 = p_to_l[q1], p_to_l[q2]
                    if l1 not in data_log and l2 not in data_log:
                        if l2 in fresh_anc and l1 not in fresh_anc:
                            active_carriers[l2] = set(active_carriers[l1])
                            fresh_anc.discard(l2)
                        elif l1 in fresh_anc and l2 not in fresh_anc:
                            active_carriers[l1] = set(active_carriers[l2])
                            fresh_anc.discard(l1)
                        else:
                            comb = active_carriers[l1] | active_carriers[l2]
                            active_carriers[l1] = set(comb)
                            active_carriers[l2] = set(comb)
                            fresh_anc.discard(l1)
                            fresh_anc.discard(l2)
                    elif l1 in data_log and l2 not in data_log:
                        fresh_anc.discard(l2)
                    elif l2 in data_log and l1 not in data_log:
                        fresh_anc.discard(l1)
            elif name == "SWAP":
                for idx in range(0, len(t_vals), 2):
                    q1, q2 = t_vals[idx], t_vals[idx + 1]
                    p_to_l[q1], p_to_l[q2] = p_to_l[q2], p_to_l[q1]
            elif name in ["CXSWAP", "CZSWAP", "ISWAP", "ISWAP_DAG"]:
                for idx in range(0, len(t_vals), 2):
                    q1, q2 = t_vals[idx], t_vals[idx + 1]
                    l1, l2 = p_to_l[q1], p_to_l[q2]
                    if l1 in data_log and l2 not in data_log:
                        fresh_anc.discard(l2)
                    elif l2 in data_log and l1 not in data_log:
                        fresh_anc.discard(l1)
                    p_to_l[q1], p_to_l[q2] = p_to_l[q2], p_to_l[q1]

    pass2(circuit, dict(phys_to_log))
    return meas_carriers


def _rewrite_circuit_internal(
    circ: Any,
    phys_to_log: Dict[int, int],
    active_gates_count: Dict[int, int]
) -> Any:
    """Rewrites a Stim circuit by tracking SWAPs virtually."""
    stim = _get_stim()
    new_circ = stim.Circuit()
    non_active_instructions = {
        "QUBIT_COORDS", "SHIFT_COORDS", "TICK",
        "DETECTOR", "OBSERVABLE_INCLUDE", "MPAD"
    }

    for inst in circ:
        if isinstance(inst, stim.CircuitRepeatBlock):
            rep_count = inst.repeat_count
            body = inst.body_copy()
            map_copy = dict(phys_to_log)
            for sub_inst in body:
                if not isinstance(sub_inst, stim.CircuitRepeatBlock):
                    if sub_inst.name in [
                        "SWAP", "CXSWAP", "CZSWAP", "ISWAP", "ISWAP_DAG"
                    ]:
                        targets = sub_inst.targets_copy()
                        for idx in range(0, len(targets), 2):
                            q1 = targets[idx].value
                            q2 = targets[idx + 1].value
                            map_copy[q1], map_copy[q2] = (
                                map_copy[q2], map_copy[q1]
                            )
            if map_copy == phys_to_log:
                new_body = _rewrite_circuit_internal(
                    body, phys_to_log, active_gates_count
                )
                new_circ.append(stim.CircuitRepeatBlock(rep_count, new_body))
            else:
                for _ in range(rep_count):
                    unrolled_body = _rewrite_circuit_internal(
                        body, phys_to_log, active_gates_count
                    )
                    new_circ += unrolled_body
            continue

        name = inst.name
        targets = inst.targets_copy()
        args = inst.gate_args_copy()

        if name == "SWAP":
            for idx in range(0, len(targets), 2):
                q1 = targets[idx].value
                q2 = targets[idx + 1].value
                phys_to_log[q1], phys_to_log[q2] = (
                    phys_to_log[q2], phys_to_log[q1]
                )
        elif name == "CXSWAP":
            cx_targets = []
            for idx in range(0, len(targets), 2):
                q1 = targets[idx].value
                q2 = targets[idx + 1].value
                l1 = phys_to_log[q1]
                l2 = phys_to_log[q2]
                cx_targets.extend([l1, l2])
                active_gates_count[l1] = active_gates_count.get(l1, 0) + 1
                active_gates_count[l2] = active_gates_count.get(l2, 0) + 1
                phys_to_log[q1], phys_to_log[q2] = (
                    phys_to_log[q2], phys_to_log[q1]
                )
            new_circ.append("CX", cx_targets, args)
        elif name == "CZSWAP":
            cz_targets = []
            for idx in range(0, len(targets), 2):
                q1 = targets[idx].value
                q2 = targets[idx + 1].value
                l1 = phys_to_log[q1]
                l2 = phys_to_log[q2]
                cz_targets.extend([l1, l2])
                active_gates_count[l1] = active_gates_count.get(l1, 0) + 1
                active_gates_count[l2] = active_gates_count.get(l2, 0) + 1
                phys_to_log[q1], phys_to_log[q2] = (
                    phys_to_log[q2], phys_to_log[q1]
                )
            new_circ.append("CZ", cz_targets, args)
        elif name in ["ISWAP", "ISWAP_DAG"]:
            cz_targets = []
            s_targets = []
            for idx in range(0, len(targets), 2):
                q1 = targets[idx].value
                q2 = targets[idx + 1].value
                l1 = phys_to_log[q1]
                l2 = phys_to_log[q2]
                cz_targets.extend([l1, l2])
                s_targets.extend([l1, l2])
                active_gates_count[l1] = active_gates_count.get(l1, 0) + 1
                active_gates_count[l2] = active_gates_count.get(l2, 0) + 1
                phys_to_log[q1], phys_to_log[q2] = (
                    phys_to_log[q2], phys_to_log[q1]
                )
            new_circ.append("CZ", cz_targets)
            s_gate = "S" if name == "ISWAP" else "S_DAG"
            new_circ.append(s_gate, s_targets)
        else:
            new_targets = []
            for t in targets:
                if t.is_qubit_target:
                    l_q = phys_to_log[t.value]
                    if name not in non_active_instructions:
                        active_gates_count[l_q] = (
                            active_gates_count.get(l_q, 0) + 1
                        )
                    if t.is_inverted_result_target:
                        new_targets.append(stim.target_inv(l_q))
                    else:
                        new_targets.append(l_q)
                elif t.is_x_target:
                    l_q = phys_to_log[t.value]
                    if name not in non_active_instructions:
                        active_gates_count[l_q] = (
                            active_gates_count.get(l_q, 0) + 1
                        )
                    new_targets.append(stim.target_x(l_q))
                elif t.is_y_target:
                    l_q = phys_to_log[t.value]
                    if name not in non_active_instructions:
                        active_gates_count[l_q] = (
                            active_gates_count.get(l_q, 0) + 1
                        )
                    new_targets.append(stim.target_y(l_q))
                elif t.is_z_target:
                    l_q = phys_to_log[t.value]
                    if name not in non_active_instructions:
                        active_gates_count[l_q] = (
                            active_gates_count.get(l_q, 0) + 1
                        )
                    new_targets.append(stim.target_z(l_q))
                else:
                    new_targets.append(t)
            new_circ.append(name, new_targets, args)
    return new_circ


def classify_qubits_thorough(
    circuit: Any,
    verbose: bool = False,
    basis_arg: Optional[str] = None,
    filepath: Optional[Union[str, Path]] = None
) -> Dict[str, Any]:
    """Classifies qubits and checks via virtual SWAP rewriting and Pauli basis tracking.

    Tracks both ancilla stabilizer evolution and data qubit initial basis
    evolution through 1-qubit and 2-qubit Clifford gates. Supports standard
    CSS codes as well as locally basis-rotated CSS codes (such as XZZX codes)
    via 2-coloring of the check-data Pauli compatibility graph.

    Returns:
        Dictionary with keys:
        - 'x_ancillas': Sorted list of X-type (or X-sector) ancilla qubit IDs.
        - 'z_ancillas': Sorted list of Z-type (or Z-sector) ancilla qubit IDs.
        - 'data_qubits': Sorted list of data qubit IDs.
        - 'routing_qubits': Sorted list of idle/routing qubit IDs.
        - 'rotated_data_qubits': Sorted list of data qubits with rotated local basis.
        - 'is_css': True if the circuit implements a (possibly locally rotated) CSS code.
        - 'is_rotated_css': True if CSS only after local data basis rotations (e.g. XZZX).
        - 'basis': Primary memory basis ('X' or 'Z').
        - 'rewritten_circuit': Circuit with SWAPs virtually eliminated.
    """
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        if filepath is None:
            filepath = circuit
        circuit = stim.Circuit.from_file(str(circuit))

    all_qubits_set: Set[int] = set()

    def collect_qubits(blk: Any) -> None:
        for inst in blk:
            if isinstance(inst, stim.CircuitRepeatBlock):
                collect_qubits(inst.body_copy())
            else:
                for t in inst.targets_copy():
                    if (
                        t.is_qubit_target or t.is_x_target
                        or t.is_y_target or t.is_z_target
                    ):
                        all_qubits_set.add(t.value)

    collect_qubits(circuit)
    max_q = max(all_qubits_set) if all_qubits_set else -1
    phys_to_log = {q: q for q in range(max_q + 1)}
    active_gates_count = {q: 0 for q in all_qubits_set}

    def cap_repeats(blk: Any) -> Any:
        out = stim.Circuit()
        for inst in blk:
            if isinstance(inst, stim.CircuitRepeatBlock):
                out.append(
                    stim.CircuitRepeatBlock(
                        min(inst.repeat_count, 2),
                        cap_repeats(inst.body_copy())
                    )
                )
            else:
                out.append(inst)
        return out

    capped_circuit = cap_repeats(circuit)
    rewritten_circ = _rewrite_circuit_internal(
        capped_circuit, phys_to_log, active_gates_count
    )
    routing_qubits = {
        q for q, count in active_gates_count.items() if count == 0
    }
    active_qubits = all_qubits_set - routing_qubits

    flat_circ = rewritten_circ.flattened()

    tot_resets: Dict[int, int] = {q: 0 for q in active_qubits}
    tot_meas: Dict[int, int] = {q: 0 for q in active_qubits}
    meas_records: List[int] = []
    used_in_detector: Set[int] = set()
    used_in_observable: Set[int] = set()
    first_detector_seen = False
    meas_before_first_det: Set[int] = set()
    r1_single_det_anc: Set[int] = set()
    multi_det_groups: List[Set[int]] = []

    for inst in flat_circ:
        name = inst.name
        targets = inst.targets_copy()
        if name in ["R", "RX", "RY", "RZ"]:
            for t in targets:
                if t.is_qubit_target and t.value in tot_resets:
                    tot_resets[t.value] += 1
        elif name in ["M", "MX", "MY", "MZ"]:
            for t in targets:
                if t.is_qubit_target:
                    if t.value in tot_meas:
                        tot_meas[t.value] += 1
                    meas_records.append(t.value)
                    if not first_detector_seen:
                        meas_before_first_det.add(t.value)
        elif name in ["MR", "MRX", "MRY", "MRZ"]:
            for t in targets:
                if t.is_qubit_target:
                    if t.value in tot_resets:
                        tot_resets[t.value] += 1
                    if t.value in tot_meas:
                        tot_meas[t.value] += 1
                    meas_records.append(t.value)
                    if not first_detector_seen:
                        meas_before_first_det.add(t.value)
        elif name == "DETECTOR":
            first_detector_seen = True
            recs = [
                meas_records[len(meas_records) + t.value]
                for t in targets
                if t.is_measurement_record_target
                and 0 <= len(meas_records) + t.value < len(meas_records)
            ]
            for q in recs:
                used_in_detector.add(q)
            if len(recs) == 1:
                r1_single_det_anc.add(recs[0])
            elif len(recs) > 1:
                multi_det_groups.append(set(recs))
        elif name == "OBSERVABLE_INCLUDE":
            for t in targets:
                if t.is_measurement_record_target:
                    idx = len(meas_records) + t.value
                    if 0 <= idx < len(meas_records):
                        used_in_observable.add(meas_records[idx])

    max_resets = max(tot_resets.values()) if tot_resets else 0
    max_meas = max(tot_meas.values()) if tot_meas else 0

    data_qubits: Set[int] = set()
    anc_candidates: Set[int] = set()

    for q in active_qubits:
        if max_resets > 1:
            is_data = (tot_resets[q] <= 1)
        elif max_meas > 1:
            is_data = (tot_meas[q] <= 1)
        else:
            if q not in used_in_detector:
                is_data = True
            else:
                is_data = (
                    (q not in meas_before_first_det)
                    or (q in used_in_observable)
                )
        if is_data:
            data_qubits.add(q)
        else:
            anc_candidates.add(q)

    # First-round stabilizer & data-basis tracking
    reset_basis: Dict[int, str] = {}
    meas_basis: Dict[int, str] = {}
    first_round_insts: List[Tuple[str, List[int]]] = []
    measured_anc: Set[int] = set()
    has_any_anc_meas = False

    for inst in flat_circ:
        name = inst.name
        t_vals = [t.value for t in inst.targets_copy() if t.is_qubit_target]
        if name in ["R", "RZ", "RX", "RY"]:
            if has_any_anc_meas and any(
                q in anc_candidates and q in measured_anc for q in t_vals
            ):
                break
            b = "X" if name == "RX" else ("Y" if name == "RY" else "Z")
            for q in t_vals:
                if q not in reset_basis:
                    reset_basis[q] = b
            first_round_insts.append((name, t_vals))
        elif name in ["M", "MZ", "MX", "MY", "MR", "MRZ", "MRX", "MRY"]:
            if any(q in data_qubits for q in t_vals) and has_any_anc_meas:
                break
            if any(q in anc_candidates and q in measured_anc for q in t_vals):
                break
            b = (
                "X" if name in ["MX", "MRX"]
                else ("Y" if name in ["MY", "MRY"] else "Z")
            )
            for q in t_vals:
                if q in anc_candidates:
                    meas_basis[q] = b
                    measured_anc.add(q)
                    has_any_anc_meas = True
                    if q not in reset_basis and name in [
                        "MR", "MRZ", "MRX", "MRY"
                    ]:
                        reset_basis[q] = b
            first_round_insts.append((name, t_vals))
        elif name in [
            "H", "H_XZ", "H_XY", "H_YZ", "S", "S_DAG",
            "SQRT_Z", "SQRT_Z_DAG", "SQRT_X", "SQRT_X_DAG",
            "SQRT_Y", "SQRT_Y_DAG",
            "CX", "ZCX", "CZ", "ZCZ", "CY", "ZCY",
            "XCX", "XCY", "XCZ", "YCX", "YCY", "YCZ"
        ]:
            first_round_insts.append((name, t_vals))

    def conj_1q(gate: str, p: str) -> str:
        if p == "I":
            return "I"
        if gate in ["H", "H_XZ"]:
            return "Z" if p == "X" else ("X" if p == "Z" else "Y")
        if gate in ["S", "S_DAG", "SQRT_Z", "SQRT_Z_DAG", "H_XY"]:
            return "Y" if p == "X" else ("X" if p == "Y" else "Z")
        if gate in ["SQRT_X", "SQRT_X_DAG", "H_YZ"]:
            return "Z" if p == "Y" else ("Y" if p == "Z" else "X")
        if gate in ["SQRT_Y", "SQRT_Y_DAG"]:
            return "Z" if p == "X" else ("X" if p == "Z" else "Y")
        return p

    mul_tab = {
        ("I", "I"): "I", ("I", "X"): "X", ("I", "Y"): "Y", ("I", "Z"): "Z",
        ("X", "I"): "X", ("X", "X"): "I", ("X", "Y"): "Z", ("X", "Z"): "Y",
        ("Y", "I"): "Y", ("Y", "X"): "Z", ("Y", "Y"): "I", ("Y", "Z"): "X",
        ("Z", "I"): "Z", ("Z", "X"): "Y", ("Z", "Y"): "X", ("Z", "Z"): "I",
    }

    def conj_2q(gate: str, p1: str, p2: str) -> Tuple[str, str]:
        if p1 == "I" and p2 == "I":
            return "I", "I"
        if gate in ["CX", "ZCX"]:
            a, b = "Z", "X"
        elif gate in ["CZ", "ZCZ"]:
            a, b = "Z", "Z"
        elif gate in ["CY", "ZCY"]:
            a, b = "Z", "Y"
        elif gate == "XCX":
            a, b = "X", "X"
        elif gate == "XCY":
            a, b = "X", "Y"
        elif gate == "XCZ":
            a, b = "X", "Z"
        elif gate == "YCX":
            a, b = "Y", "X"
        elif gate == "YCY":
            a, b = "Y", "Y"
        elif gate == "YCZ":
            a, b = "Y", "Z"
        else:
            return p1, p2
        anti1 = (p1 != "I" and p1 != a)
        anti2 = (p2 != "I" and p2 != b)
        np1 = mul_tab[(p1, a)] if anti2 else p1
        np2 = mul_tab[(p2, b)] if anti1 else p2
        return np1, np2

    # Simultaneously evolve all ancilla stabilizers and data initial bases
    # using bitmasks over sorted anc_list for O(1) gate updates.
    anc_list = sorted(list(anc_candidates))
    anc_to_bit = {anc: (1 << i) for i, anc in enumerate(anc_list)}
    all_anc_mask = (1 << len(anc_list)) - 1
    active_anc_mask = 0
    qx: Dict[int, int] = {q: 0 for q in all_qubits_set}
    qz: Dict[int, int] = {q: 0 for q in all_qubits_set}

    for anc in anc_list:
        b = reset_basis.get(anc, meas_basis.get(anc, "Z"))
        bit = anc_to_bit[anc]
        if b in ["X", "Y"]:
            qx[anc] |= bit
        if b in ["Z", "Y"]:
            qz[anc] |= bit
        if anc not in reset_basis:
            active_anc_mask |= bit

    data_init_basis: Dict[int, str] = {q: "Z" for q in data_qubits}

    for name, t_vals in first_round_insts:
        if name in ["R", "RZ", "RX", "RY"]:
            b = "X" if name == "RX" else ("Y" if name == "RY" else "Z")
            for q in t_vals:
                if q in data_qubits:
                    data_init_basis[q] = b
                if q in anc_to_bit:
                    bit = anc_to_bit[q]
                    if not (active_anc_mask & bit):
                        active_anc_mask |= bit
                        if b in ["X", "Y"]:
                            qx[q] |= bit
                        else:
                            qx[q] &= ~bit
                        if b in ["Z", "Y"]:
                            qz[q] |= bit
                        else:
                            qz[q] &= ~bit
        elif name in ["M", "MZ", "MX", "MY", "MR", "MRZ", "MRX", "MRY"]:
            for q in t_vals:
                if q in anc_to_bit:
                    active_anc_mask &= ~anc_to_bit[q]
        elif name in [
            "H", "H_XZ", "H_XY", "H_YZ", "S", "S_DAG",
            "SQRT_Z", "SQRT_Z_DAG", "SQRT_X", "SQRT_X_DAG",
            "SQRT_Y", "SQRT_Y_DAG"
        ]:
            for q in t_vals:
                if q in data_qubits:
                    data_init_basis[q] = conj_1q(name, data_init_basis[q])
                if q not in qx:
                    continue
                mask = active_anc_mask if q in anc_candidates else all_anc_mask
                x_m = qx[q] & mask
                z_m = qz[q] & mask
                if name in ["H", "H_XZ", "SQRT_Y", "SQRT_Y_DAG"]:
                    nx_m, nz_m = z_m, x_m
                elif name in [
                    "S", "S_DAG", "SQRT_Z", "SQRT_Z_DAG", "H_XY"
                ]:
                    nx_m, nz_m = x_m, x_m ^ z_m
                elif name in ["SQRT_X", "SQRT_X_DAG", "H_YZ"]:
                    nx_m, nz_m = x_m ^ z_m, z_m
                else:
                    nx_m, nz_m = x_m, z_m
                qx[q] = (qx[q] & ~mask) | nx_m
                qz[q] = (qz[q] & ~mask) | nz_m
        elif name in [
            "CX", "ZCX", "CZ", "ZCZ", "CY", "ZCY",
            "XCX", "XCY", "XCZ", "YCX", "YCY", "YCZ"
        ]:
            if name in ["CX", "ZCX"]:
                ax, az, bx, bz = 0, all_anc_mask, all_anc_mask, 0
            elif name in ["CZ", "ZCZ"]:
                ax, az, bx, bz = 0, all_anc_mask, 0, all_anc_mask
            elif name in ["CY", "ZCY"]:
                ax, az, bx, bz = 0, all_anc_mask, all_anc_mask, all_anc_mask
            elif name == "XCX":
                ax, az, bx, bz = all_anc_mask, 0, all_anc_mask, 0
            elif name == "XCY":
                ax, az, bx, bz = all_anc_mask, 0, all_anc_mask, all_anc_mask
            elif name == "XCZ":
                ax, az, bx, bz = all_anc_mask, 0, 0, all_anc_mask
            elif name == "YCX":
                ax, az, bx, bz = all_anc_mask, all_anc_mask, all_anc_mask, 0
            elif name == "YCY":
                ax, az, bx, bz = (
                    all_anc_mask, all_anc_mask, all_anc_mask, all_anc_mask
                )
            elif name == "YCZ":
                ax, az, bx, bz = all_anc_mask, all_anc_mask, 0, all_anc_mask
            else:
                continue

            for idx in range(0, len(t_vals), 2):
                q1, q2 = t_vals[idx], t_vals[idx + 1]
                if q1 not in qx or q2 not in qx:
                    continue
                mask = active_anc_mask
                if q1 in anc_to_bit and not (active_anc_mask & anc_to_bit[q1]):
                    mask &= ~anc_to_bit[q1]
                if q2 in anc_to_bit and not (active_anc_mask & anc_to_bit[q2]):
                    mask &= ~anc_to_bit[q2]
                x1 = qx[q1] & mask
                z1 = qz[q1] & mask
                x2 = qx[q2] & mask
                z2 = qz[q2] & mask
                anti1 = (x1 & (az & mask)) ^ (z1 & (ax & mask))
                anti2 = (x2 & (bz & mask)) ^ (z2 & (bx & mask))
                nx1 = x1 ^ (anti2 & ax)
                nz1 = z1 ^ (anti2 & az)
                nx2 = x2 ^ (anti1 & bx)
                nz2 = z2 ^ (anti1 & bz)
                qx[q1] = (qx[q1] & ~mask) | nx1
                qz[q1] = (qz[q1] & ~mask) | nz1
                qx[q2] = (qx[q2] & ~mask) | nx2
                qz[q2] = (qz[q2] & ~mask) | nz2

    anc_data_paulis: Dict[int, Dict[int, str]] = {
        anc: {} for anc in anc_list
    }
    for q in data_qubits:
        xm = qx[q]
        zm = qz[q]
        active_m = xm | zm
        while active_m:
            lsb = active_m & -active_m
            anc_idx = lsb.bit_length() - 1
            anc = anc_list[anc_idx]
            has_x = bool(xm & lsb)
            has_z = bool(zm & lsb)
            p_str = (
                "Y" if (has_x and has_z)
                else ("X" if has_x else "Z")
            )
            anc_data_paulis[anc][q] = p_str
            active_m ^= lsb

    active_ancillas = [a for a in anc_list if len(anc_data_paulis[a]) > 0]
    idle_ancillas = [a for a in anc_list if len(anc_data_paulis[a]) == 0]
    routing_qubits.update(idle_ancillas)

    # Check commutation among active ancillas
    checks_commute = True
    for i in range(len(active_ancillas)):
        a1 = active_ancillas[i]
        s1 = anc_data_paulis[a1]
        for j in range(i + 1, len(active_ancillas)):
            a2 = active_ancillas[j]
            s2 = anc_data_paulis[a2]
            anti = 0
            for q, p1 in s1.items():
                p2 = s2.get(q)
                if p2 is not None and p1 != p2:
                    anti ^= 1
            if anti != 0:
                checks_commute = False
                break
        if not checks_commute:
            break

    # Build check-relation graph on active_ancillas via shared data qubits
    adj: Dict[int, Dict[int, int]] = {a: {} for a in active_ancillas}
    bipartite_paulis = True
    for q in data_qubits:
        touching = [
            (a, anc_data_paulis[a][q])
            for a in active_ancillas
            if q in anc_data_paulis[a]
        ]
        if not touching:
            continue
        pauli_set = {p for _, p in touching}
        if len(pauli_set) > 2:
            bipartite_paulis = False
            break
        a0, p0 = touching[0]
        for ak, pk in touching[1:]:
            rel = 0 if pk == p0 else 1
            if ak in adj[a0] and adj[a0][ak] != rel:
                bipartite_paulis = False
                break
            adj[a0][ak] = rel
            adj[ak][a0] = rel
        if not bipartite_paulis:
            break

    # 2-color connected components of active_ancillas
    color: Dict[int, int] = {}
    components: List[Tuple[List[int], List[int]]] = []
    if bipartite_paulis:
        for a in active_ancillas:
            if a in color:
                continue
            c0: List[int] = []
            c1: List[int] = []
            color[a] = 0
            queue = [a]
            idx_q = 0
            while idx_q < len(queue):
                u = queue[idx_q]
                idx_q += 1
                if color[u] == 0:
                    c0.append(u)
                else:
                    c1.append(u)
                for v, rel in adj[u].items():
                    expected = color[u] ^ rel
                    if v not in color:
                        color[v] = expected
                        queue.append(v)
                    elif color[v] != expected:
                        bipartite_paulis = False
                        break
                if not bipartite_paulis:
                    break
            if not bipartite_paulis:
                break
            components.append((c0, c1))

    # Determine if all checks are already pure-X or pure-Z (standard CSS)
    all_pure_xz = True
    for a in active_ancillas:
        ps = set(anc_data_paulis[a].values())
        if not (ps == {"X"} or ps == {"Z"}):
            all_pure_xz = False
            break

    # Detect color code (pure X and pure Z checks with identical support)
    native_basis = detect_basis(
        circuit, filepath=filepath, basis_arg=basis_arg
    )
    final_det_anc: Set[int] = set()
    for grp in multi_det_groups:
        if grp & data_qubits:
            final_det_anc.update(grp & set(active_ancillas))

    x_ancillas: Set[int] = set()
    z_ancillas: Set[int] = set()
    rotated_data_qubits: Set[int] = set()
    is_css = False
    is_rotated_css = False

    if checks_commute and bipartite_paulis and active_ancillas:
        if all_pure_xz:
            for a in active_ancillas:
                ps = set(anc_data_paulis[a].values())
                if ps == {"X"}:
                    x_ancillas.add(a)
                else:
                    z_ancillas.add(a)
            x_supps = {
                frozenset(anc_data_paulis[a].keys()) for a in x_ancillas
            }
            z_supps = {
                frozenset(anc_data_paulis[a].keys()) for a in z_ancillas
            }
            is_color_code = (
                len(x_supps) > 0 and x_supps == z_supps
                and len(x_ancillas) == len(z_ancillas)
            )
            is_css = not is_color_code
            is_rotated_css = False
        else:
            is_css = True
            is_rotated_css = True
            for c0, c1 in components:
                # Score which sector is the primary memory basis
                def score_primary(sector: List[int]) -> Tuple[int, int, int]:
                    r1_cnt = sum(1 for a in sector if a in r1_single_det_anc)
                    fin_cnt = sum(1 for a in sector if a in final_det_anc)
                    init_cnt = sum(
                        1 for a in sector
                        if all(
                            anc_data_paulis[a][q] == data_init_basis.get(q, "Z")
                            for q in anc_data_paulis[a]
                        )
                    )
                    return (r1_cnt, fin_cnt, init_cnt)

                sc0 = score_primary(c0)
                sc1 = score_primary(c1)
                if sc0 >= sc1:
                    primary_sec, minority_sec = c0, c1
                else:
                    primary_sec, minority_sec = c1, c0

                if native_basis == "Z":
                    z_ancillas.update(primary_sec)
                    x_ancillas.update(minority_sec)
                else:
                    x_ancillas.update(primary_sec)
                    z_ancillas.update(minority_sec)

            # Determine rotated data qubits (where check Pauli != sector label)
            for q in data_qubits:
                for a in active_ancillas:
                    if q in anc_data_paulis[a]:
                        sector_label = "X" if a in x_ancillas else "Z"
                        if anc_data_paulis[a][q] != sector_label:
                            rotated_data_qubits.add(q)
                        break
    else:
        # Non-CSS fallback classification by majority Pauli weight
        is_css = False
        is_rotated_css = False
        for a in active_ancillas:
            s = anc_data_paulis[a]
            nx = sum(1 for p in s.values() if p == "X")
            nz = sum(1 for p in s.values() if p == "Z")
            if nz > nx:
                z_ancillas.add(a)
            else:
                x_ancillas.add(a)

    return {
        "x_ancillas": sorted(list(x_ancillas)),
        "z_ancillas": sorted(list(z_ancillas)),
        "data_qubits": sorted(list(data_qubits)),
        "routing_qubits": sorted(list(routing_qubits)),
        "rotated_data_qubits": sorted(list(rotated_data_qubits)),
        "is_css": is_css,
        "is_rotated_css": is_rotated_css,
        "basis": native_basis,
        "rewritten_circuit": rewritten_circ,
    }


def verify_qubit_metadata(
    circuit: Any,
    verbose: bool = False,
    basis_arg: Optional[str] = None,
    filepath: Optional[Union[str, Path]] = None
) -> Dict[str, Any]:
    """Verifies qubit classification and returns basis tracking metadata."""
    return classify_qubits_thorough(
        circuit, verbose=verbose, basis_arg=basis_arg, filepath=filepath
    )


def strip_minority_detectors(
    circuit: Any,
    basis: str,
    verbose: bool = False,
    thorough: bool = True,
    thorough_res: Optional[Dict[str, Any]] = None
) -> Tuple[Any, int, int]:
    """Strips detectors corresponding to the minority basis from a Stim circuit.

    Always uses ancilla and data basis tracking (and dynamic carrier tracking
    across SWAPs and syndrome handoffs) to identify each detector's basis.

    Args:
        circuit: A stim.Circuit object or path to a .stim file.
        basis: Primary memory basis ('X' or 'Z'). Detectors of the opposite
            basis are removed.
        verbose: Whether to print warnings/statistics.
        thorough: Kept for API compatibility (basis tracking is always enabled).
        thorough_res: Optional precomputed result from classify_qubits_thorough.

    Returns:
        Tuple of (simplified_circuit, stripped_count, kept_count).
    """
    stim = _get_stim()
    if isinstance(circuit, (str, Path)):
        circuit = stim.Circuit.from_file(str(circuit))

    t_res = (
        thorough_res
        if thorough_res is not None
        else classify_qubits_thorough(circuit, verbose=False, basis_arg=basis)
    )
    x_anc = set(t_res["x_ancillas"])
    z_anc = set(t_res["z_ancillas"])
    minority_basis = "Z" if basis.upper() == "X" else "X"

    meas_carriers = track_circuit_carriers(circuit)

    # Coordinate metadata for consistency check on standard (non-rotated) CSS
    qubit_coords: Dict[int, Tuple[float, ...]] = {}

    def get_coords(blk: Any) -> None:
        for inst in blk:
            if isinstance(inst, stim.CircuitRepeatBlock):
                get_coords(inst.body_copy())
            elif inst.name == "QUBIT_COORDS":
                args = inst.gate_args_copy()
                for t in inst.targets_copy():
                    if t.is_qubit_target:
                        qubit_coords[t.value] = tuple(float(x) for x in args)

    get_coords(circuit)
    x_coords = {
        qubit_coords[q][:2]
        for q in x_anc
        if q in qubit_coords and len(qubit_coords[q]) >= 2
    }
    z_coords = {
        qubit_coords[q][:2]
        for q in z_anc
        if q in qubit_coords and len(qubit_coords[q]) >= 2
    }
    layout_well_defined = (
        not t_res.get("is_rotated_css", False)
        and len(x_coords) > 0
        and len(z_coords) > 0
        and x_coords.isdisjoint(z_coords)
    )

    meas_history: List[int] = []
    stripped_count = 0
    kept_count = 0
    contradictions_c_count = 0
    floquet_suspected = False

    def process_block(block: Any) -> Any:
        nonlocal stripped_count, kept_count
        nonlocal contradictions_c_count, floquet_suspected
        new_block = stim.Circuit()
        for inst in block:
            if isinstance(inst, stim.CircuitRepeatBlock):
                rep_count = inst.repeat_count
                if rep_count <= 0:
                    continue
                body_copy = inst.body_copy()
                meas_before = len(meas_history)
                strip_before = stripped_count
                kept_before = kept_count
                new_body = process_block(body_copy)
                meas_after = len(meas_history)
                meas_per_rep = meas_after - meas_before
                if rep_count > 1:
                    stripped_count += (stripped_count - strip_before) * (
                        rep_count - 1
                    )
                    kept_count += (kept_count - kept_before) * (rep_count - 1)
                    if meas_per_rep > 0:
                        meas_history.extend(
                            meas_history[meas_before:meas_after]
                            * (rep_count - 1)
                        )
                new_block.append(stim.CircuitRepeatBlock(rep_count, new_body))
            elif inst.name in [
                "M", "MX", "MY", "MZ", "MR", "MRX", "MRY", "MRZ"
            ]:
                if inst.name in ["MY", "MRY"]:
                    floquet_suspected = True
                for t in inst.targets_copy():
                    if t.is_qubit_target:
                        meas_history.append(t.value)
                new_block.append(inst)
            elif inst.name == "DETECTOR":
                cur_meas = len(meas_history)
                rec_targets = [
                    t for t in inst.targets_copy()
                    if t.is_measurement_record_target
                ]
                if not rec_targets:
                    stripped_count += 1
                    continue

                qubits = [
                    meas_history[cur_meas + t.value]
                    for t in rec_targets
                    if 0 <= cur_meas + t.value < cur_meas
                ]
                is_x_qubit = any(q in x_anc for q in qubits)
                is_z_qubit = any(q in z_anc for q in qubits)
                class_qubit = (
                    "Both" if (is_x_qubit and is_z_qubit)
                    else "X" if is_x_qubit
                    else "Z" if is_z_qubit
                    else "Unknown"
                )

                args = inst.gate_args_copy()
                class_b = "Unknown"
                if layout_well_defined and args and len(args) >= 2:
                    det_xy = (float(args[0]), float(args[1]))
                    if det_xy in x_coords:
                        class_b = "X"
                    elif det_xy in z_coords:
                        class_b = "Z"

                class_c = "Unknown"
                if meas_carriers:
                    carriers = [
                        c
                        for t in rec_targets
                        if 0 <= cur_meas + t.value < len(meas_carriers)
                        for c in meas_carriers[cur_meas + t.value]
                    ]
                    is_x_c = any(c in x_anc for c in carriers)
                    is_z_c = any(c in z_anc for c in carriers)
                    class_c = (
                        "Both" if (is_x_c and is_z_c)
                        else "X" if is_x_c
                        else "Z" if is_z_c
                        else "Unknown"
                    )

                is_minority = False
                if class_c != "Unknown":
                    if class_c == "Both":
                        contradictions_c_count += 1
                        floquet_suspected = True
                    elif class_b != "Unknown" and class_c != class_b:
                        contradictions_c_count += 1
                    elif len(rec_targets) == 1 and class_c == minority_basis:
                        contradictions_c_count += 1
                    is_minority = (class_c == minority_basis)
                elif layout_well_defined and class_b != "Unknown":
                    is_minority = (class_b == minority_basis)
                elif len(rec_targets) == 1:
                    is_minority = False
                else:
                    if class_qubit == "Both":
                        floquet_suspected = True
                    is_minority = (class_qubit == minority_basis)

                if is_minority:
                    stripped_count += 1
                else:
                    kept_count += 1
                    new_block.append(inst)
            else:
                new_block.append(inst)
        return new_block

    simp_circuit = process_block(circuit)

    if contradictions_c_count > 0 and verbose:
        sys.stderr.write(
            f"# Warning: {contradictions_c_count} detector contradiction(s) "
            "detected during dynamic carrier tracking.\n"
        )
    if floquet_suspected and verbose:
        sys.stderr.write(
            "# Warning: Floquet or dynamic measurement structure suspected; "
            "simplified circuit distance may be invalid.\n"
        )

    return simp_circuit, stripped_count, kept_count


def compute_classical_distance(
    H: Any,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmin: int = 1,
    wmax: int = 0,
    smax: Optional[int] = None,
    start: Optional[Union[int, str, Sequence[int]]] = None,
    cbeg: Optional[int] = None,
    cend: Optional[int] = None,
    noscan: int = 0,
    dW: int = -1,
    maxC: int = 0,
    finC: Optional[str] = None,
    outC: Optional[str] = None,
    do_cws: bool = False,
    return_info: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
    solver: str = "dist_m4ri",
    codedistance_method: str = "QDistEvol",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0,
    verbose: bool = False,
    nothrottle: bool = False,
    chunk_size: int = 0,
    ksub: int = 0,
    kwin: int = 0,
    win_mode: int = 0,
    min_hits: Optional[int] = None,
    cov_cws: int = 100,
    refresh: int = 0,
    trust_start: bool = False,
    **kwargs
) -> Any:
    """
    Computes the minimum distance of a classical linear code given parity check matrix H.

    Args:
        H: Parity check matrix (numpy array, scipy sparse matrix, or file path).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance.
        d_max / dmax: Known upper bound on distance.
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        smax: Maximum syndrome weight for CC confinement profile.
        start: Expert option: list of CC start columns (int, comma-separated string "a,b,c", or
            sequence of ints; None or -1: all columns).  CC clusters are grown only from the listed
            columns, and they are not limited to larger column indices.  WARNING: the lower bound
            (and an exact result) is valid only if a code symmetry maps every minimum-weight codeword
            to one containing a listed column (e.g., one column per block of a quasi-cyclic code);
            the upper bound is always valid.  Results are cached under the separate key
            '<key>:start=a,b,c'; only dmax and codewords update the main cache record.
        trust_start: Expert option: accept the start-list results as valid, i.e., also merge dmin and
            rw_steps into the main cache record, and copy an existing exact start-list record to the
            main record even if no calculation is needed (default: False).
        cbeg / cend / noscan: Disabled in the Python interface: accepted for backward compatibility,
            but ignored with a warning to stderr (use the dist_m4ri binary directly).
        ksub: Experimental, should not be used (a warning is printed to stderr for ksub > 0 with method=1 or 3):
            subspace dimension sampled from ker(H) in RW (default: 0, full-matrix RW).
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        finC: Input file with initial codewords.
        outC: Output file to save codewords (NZLIST format).
        do_cws: Whether to return extracted codewords.
        return_info: If True, return (dist, d_info) or (dist, d_info, cws) where d_info is [dmin, dmax, num_rw].
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug bitmap: the bits in BIN_DBG_MASK for the dist_m4ri binary, and the PY_DBG_* bits.
        verbose: Verbose reporting flag.

    Returns:
        dist or (dist, cws) if do_cws is True (or (dist, d_info) / (dist, d_info, cws) if return_info=True)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max
    start_list = _prepare_start_option(start, trust_start, noscan, cbeg, cend, solver)
    _warn_experimental_ksub(ksub, method, solver)

    finC = check_finc_outc(finC, outC, verbose=verbose)

    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None
    plain_key = None
    plain_entry = None

    if solver == "codedistance":
        if do_cws or outC:
            raise ValueError("Codeword extraction is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
        import numpy as np
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        H_mat = H.toarray() if hasattr(H, 'toarray') else (
            np.asarray(H, dtype=np.int8) if isinstance(H, (np.ndarray, list)) else None
        )
        res = codedistance.codeDistance(
            H_mat, None, tB=1, method=codedistance_method, params=params,
            seed=seed if seed != 0 else None
        )
        d = res.get("d", -1)
        if return_info:
            d_info = format_bounds_list(d, d, 0) if d > 0 else [0, 0, 0]
            return d, d_info
        return d

    # Solver is native multithreaded dist_m4ri (supports bounds caching and cumulative RW steps)
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            h_state = get_sparse_array_state(H)
            plain_key = f"classical:{h_state}"
            code_key, cached_entry, plain_entry = _resolve_start_cache(
                plain_key, start_list, lambda e: _single_entry_is_exact(e, bool(do_cws or outC)), trust_start,
                lambda e: _merge_start_into_plain(plain_key, e, True), eff_cache_file, verbose
            )
            eff_dmin, eff_dmax = _seed_bounds_from_entry(eff_dmin, eff_dmax, plain_entry)
            if cached_entry is not None:
                # If exact distance is already proven and not asking for more codewords
                if cached_entry.get("dmin", 0) > 0 and cached_entry.get("dmin") == cached_entry.get("dmax"):
                    if not (do_cws or outC) or (cached_entry.get("cws") and len(cached_entry["cws"]) > 0):
                        d_info = format_bounds_list(
                            cached_entry.get("dmin", 0), cached_entry.get("dmax", 0), cached_entry.get("rw_steps", 0)
                        )
                        if verbose:
                            print(f"[dist_m4ri] Cache retrieval: SUCCESS "
                                  f"(found cached exact distance for '{code_key}')")
                            print(f"[dist_m4ri] Cached result: dist={cached_entry['dist']}, "
                                  f"bounds={format_bounds_str(d_info)}")
                        elif debug & PY_DBG_CACHE:
                            print("[dist_m4ri] Cache hit for classical distance (exact distance known)!")
                        cws_res = cached_entry.get("cws", [])
                        if outC and cws_res:
                            _write_nzlist_file(outC, cws_res)
                        if return_info:
                            return (
                                cached_entry["dist"], d_info, cws_res
                            ) if do_cws else (
                                cached_entry["dist"], d_info
                            )
                        return (cached_entry["dist"], cws_res) if do_cws else cached_entry["dist"]
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: PARTIAL (cached bounds: "
                          f"dmin={cached_entry.get('dmin', 0)}, dmax={cached_entry.get('dmax', 0)}, "
                          f"rw_steps={cached_entry.get('rw_steps', 0)}; continuing search)")
                # Use existing cached bounds to accelerate subsequent runs
                if eff_dmax == 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = cached_entry["dmax"]
                elif eff_dmax > 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = min(eff_dmax, cached_entry["dmax"])
                if eff_dmin <= 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = cached_entry["dmin"]
                elif eff_dmin > 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = max(eff_dmin, cached_entry["dmin"])
            else:
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: MISS (no entry for '{code_key}')")
        except Exception:
            code_key = None
            cached_entry = None
            plain_key = None
            plain_entry = None
    else:
        if verbose:
            print("[dist_m4ri] Cache retrieval: DISABLED (cache is turned off)")

    temp_files = []
    try:
        if isinstance(H, (str, Path)) and os.path.exists(str(H)):
            file_H = str(H)
        else:
            file_H = _matrix_to_file(H, extension="_H.mtx")
            temp_files.append(file_H)

        outC_file = None
        if do_cws or outC:
            outC_file = create_unique_file(extension="_cws.nz")
            temp_files.append(outC_file)

        dmin_res, dmax_res, rw_steps = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            finH=file_H,
            finC=finC,
            classical=1,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmin=wmin,
            wmax=wmax,
            smax=smax,
            start=start_list,
            warn_start=False,
            warn_ksub=False,
            dexp=d_exp,
            steps=num_steps,
            threads=threads,
            timeout=timeout,
            dW=dW,
            maxC=maxC,
            outC=outC_file,
            seed=seed,
            debug=debug,
            nothrottle=nothrottle,
            chunk_size=chunk_size,
            ksub=ksub,
            kwin=kwin if kwin > 0 else int(kwargs.get("win", 0) or 0),
            win_mode=win_mode,
            min_hits=min_hits,
            cov_cws=cov_cws,
            refresh=refresh,
            verbose=verbose
        )

        dist = dmin_res if (dmin_res == dmax_res or dmax_res == 0) else dmax_res
        cws = []
        if (do_cws or outC) and outC_file and os.path.exists(outC_file):
            cws = read_sparse_vectors(outC_file)
            cws.sort(key=len)
            if outC:
                _write_nzlist_file(outC, cws)

        d_info = format_bounds_list(dmin_res, dmax_res, rw_steps)

        if _use_distance_cache and code_key is not None:
            prev_steps = cached_entry.get("rw_steps", 0) if cached_entry else 0
            prev_dmax = cached_entry.get("dmax", 0) if cached_entry else 0
            prev_dmin = cached_entry.get("dmin", 0) if cached_entry else 0
            prev_cws = list(cached_entry.get("cws", [])) if cached_entry else []

            total_rw_steps = prev_steps + rw_steps
            best_dmax = (
                min(prev_dmax, dmax_res) if (prev_dmax > 0 and dmax_res > 0)
                else (dmax_res if dmax_res > 0 else prev_dmax)
            )
            best_dmin = max(prev_dmin, dmin_res)

            combined_cws = prev_cws
            if cws:
                existing_set = {tuple(cw) for cw in combined_cws}
                for cw in cws:
                    if tuple(cw) not in existing_set:
                        combined_cws.append(cw)
                        existing_set.add(tuple(cw))
                combined_cws.sort(key=len)

            d_info = format_bounds_list(best_dmin, best_dmax, total_rw_steps)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "d_info": d_info,
                "cws": combined_cws
            }
            if plain_key is not None and code_key != plain_key:
                # start-list run: merge the valid part (or all, with trust_start) into the main record
                _merge_start_into_plain(plain_key, _distance_cache[code_key], trust_start, {"": rw_steps})
            if eff_cache_file:
                save_distance_cache(eff_cache_file)

            if return_info:
                return (dist, d_info, combined_cws) if do_cws else (dist, d_info)
            return (dist, combined_cws) if do_cws else dist

        if return_info:
            return (dist, d_info, cws) if do_cws else (dist, d_info)
        return (dist, cws) if do_cws else dist

    finally:
        _remove_temp_files(temp_files, debug)


def compute_quantum_distance(
    H: Any,
    G: Optional[Any] = None,
    L: Optional[Any] = None,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmin: int = 1,
    wmax: int = 0,
    smax: Optional[int] = None,
    start: Optional[Union[int, str, Sequence[int]]] = None,
    cbeg: Optional[int] = None,
    cend: Optional[int] = None,
    noscan: int = 0,
    dW: int = -1,
    maxC: int = 0,
    finC: Optional[str] = None,
    outC: Optional[str] = None,
    do_cws: bool = False,
    return_info: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
    solver: str = "dist_m4ri",
    codedistance_method: str = "QDistEvol",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0,
    verbose: bool = False,
    nothrottle: bool = False,
    chunk_size: int = 0,
    ksub: int = 0,
    kwin: int = 0,
    win_mode: int = 0,
    min_hits: Optional[int] = None,
    cov_cws: int = 100,
    refresh: int = 0,
    trust_start: bool = False,
    **kwargs
) -> Any:
    """
    Computes the minimum distance of a single-sided quantum code given parity check matrix H
    and degeneracy generator G (or logical operator matrix L).

    Args:
        H: Parity check matrix (numpy array, scipy sparse matrix, or file path).
        G: Degeneracy generator matrix (numpy array, scipy sparse matrix, or file path).
        L: Logical operator matrix (numpy array, scipy sparse matrix, or file path).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance.
        d_max / dmax: Known upper bound on distance.
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        smax: Maximum syndrome weight for CC confinement profile.
        start: Expert option: list of CC start columns (int, comma-separated string "a,b,c", or
            sequence of ints; None or -1: all columns).  CC clusters are grown only from the listed
            columns, and they are not limited to larger column indices.  WARNING: the lower bound
            (and an exact result) is valid only if a code symmetry maps every minimum-weight codeword
            to one containing a listed column (e.g., one column per block of a quasi-cyclic code);
            the upper bound is always valid.  Results are cached under the separate key
            '<key>:start=a,b,c'; only dmax and codewords update the main cache record.
        trust_start: Expert option: accept the start-list results as valid, i.e., also merge dmin and
            rw_steps into the main cache record, and copy an existing exact start-list record to the
            main record even if no calculation is needed (default: False).
        cbeg / cend / noscan: Disabled in the Python interface: accepted for backward compatibility,
            but ignored with a warning to stderr (use the dist_m4ri binary directly).
        ksub: Experimental, should not be used (a warning is printed to stderr for ksub > 0 with method=1 or 3):
            subspace dimension sampled from ker(H) in RW (default: 0, full-matrix RW).
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        finC: Input file with initial codewords.
        outC: Output file to save codewords (NZLIST format).
        do_cws: Whether to return extracted codewords.
        return_info: If True, return (dist, d_info) or (dist, d_info, cws) where d_info is [dmin, dmax, num_rw].
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug bitmap: the bits in BIN_DBG_MASK for the dist_m4ri binary, and the PY_DBG_* bits.
        verbose: Verbose reporting flag.

    Returns:
        dist or (dist, cws) if do_cws is True (or (dist, d_info) / (dist, d_info, cws) if return_info=True)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max
    start_list = _prepare_start_option(start, trust_start, noscan, cbeg, cend, solver)
    _warn_experimental_ksub(ksub, method, solver)

    if G is None and L is None:
        raise ValueError(
            "Either G (dual generator matrix) or L (logical operator matrix) "
            "must be specified for quantum distance."
        )

    finC = check_finc_outc(finC, outC, verbose=verbose)

    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None
    plain_key = None
    plain_entry = None

    if solver == "codedistance":
        if do_cws or outC:
            raise ValueError("Codeword extraction is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
        import numpy as np
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        H_mat = H.toarray() if hasattr(H, 'toarray') else (
            np.asarray(H, dtype=np.int8) if isinstance(H, (np.ndarray, list)) else None
        )
        if G is None and L is not None:
            raise ValueError(
                "The 'codedistance' solver requires the stabilizer generator matrix G, "
                "not the logical operator matrix L. Use solver='dist_m4ri' with L."
            )
        dual_mat = G if G is not None else L
        dual_arr = dual_mat.toarray() if hasattr(dual_mat, 'toarray') else (
            np.asarray(dual_mat, dtype=np.int8) if isinstance(dual_mat, (np.ndarray, list)) else None
        )

        res = codedistance.codeDistance(
            H_mat, dual_arr, tB=1, method=codedistance_method, params=params,
            seed=seed if seed != 0 else None
        )
        d = res.get("d", -1)
        if return_info:
            d_info = format_bounds_list(d, d, 0) if d > 0 else [0, 0, 0]
            return d, d_info
        return d

    # Solver is native multithreaded dist_m4ri
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            h_state = get_sparse_array_state(H)
            if G is not None:
                g_state = get_sparse_array_state(G)
                plain_key = f"quantum:H={h_state}:G={g_state}"
            else:
                l_state = get_sparse_array_state(L)
                plain_key = f"quantum:H={h_state}:L={l_state}"
            code_key, cached_entry, plain_entry = _resolve_start_cache(
                plain_key, start_list, lambda e: _single_entry_is_exact(e, bool(do_cws or outC)), trust_start,
                lambda e: _merge_start_into_plain(plain_key, e, True), eff_cache_file, verbose
            )
            eff_dmin, eff_dmax = _seed_bounds_from_entry(eff_dmin, eff_dmax, plain_entry)
            if cached_entry is not None:
                if cached_entry.get("dmin", 0) > 0 and cached_entry.get("dmin") == cached_entry.get("dmax"):
                    if not (do_cws or outC) or (cached_entry.get("cws") and len(cached_entry["cws"]) > 0):
                        d_info = format_bounds_list(
                            cached_entry.get("dmin", 0), cached_entry.get("dmax", 0), cached_entry.get("rw_steps", 0)
                        )
                        if verbose:
                            print(f"[dist_m4ri] Cache retrieval: SUCCESS "
                                  f"(found cached exact distance for '{code_key}')")
                            print(f"[dist_m4ri] Cached result: dist={cached_entry['dist']}, "
                                  f"bounds={format_bounds_str(d_info)}")
                        elif debug & PY_DBG_CACHE:
                            print("[dist_m4ri] Cache hit for quantum distance (exact distance known)!")
                        cws_res = cached_entry.get("cws", [])
                        if outC and cws_res:
                            _write_nzlist_file(outC, cws_res)
                        if return_info:
                            return (
                                cached_entry["dist"], d_info, cws_res
                            ) if do_cws else (
                                cached_entry["dist"], d_info
                            )
                        return (cached_entry["dist"], cws_res) if do_cws else cached_entry["dist"]
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: PARTIAL (cached bounds: "
                          f"dmin={cached_entry.get('dmin', 0)}, dmax={cached_entry.get('dmax', 0)}, "
                          f"rw_steps={cached_entry.get('rw_steps', 0)}; continuing search)")
                if eff_dmax == 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = cached_entry["dmax"]
                elif eff_dmax > 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = min(eff_dmax, cached_entry["dmax"])
                if eff_dmin <= 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = cached_entry["dmin"]
                elif eff_dmin > 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = max(eff_dmin, cached_entry["dmin"])
            else:
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: MISS (no entry for '{code_key}')")
        except Exception:
            code_key = None
            cached_entry = None
            plain_key = None
            plain_entry = None
    else:
        if verbose:
            print("[dist_m4ri] Cache retrieval: DISABLED (cache is turned off)")

    temp_files = []
    try:
        if isinstance(H, (str, Path)) and os.path.exists(str(H)):
            file_H = str(H)
        else:
            file_H = _matrix_to_file(H, extension="_H.mtx")
            temp_files.append(file_H)

        file_G = None
        if G is not None:
            if isinstance(G, (str, Path)) and os.path.exists(str(G)):
                file_G = str(G)
            else:
                file_G = _matrix_to_file(G, extension="_G.mtx")
                temp_files.append(file_G)

        file_L = None
        if L is not None:
            if isinstance(L, (str, Path)) and os.path.exists(str(L)):
                file_L = str(L)
            else:
                file_L = _matrix_to_file(L, extension="_L.mtx")
                temp_files.append(file_L)

        outC_file = None
        if do_cws or outC:
            outC_file = create_unique_file(extension="_cws.nz")
            temp_files.append(outC_file)

        dmin_res, dmax_res, rw_steps = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            finH=file_H,
            finG=file_G,
            finL=file_L,
            finC=finC,
            classical=0,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmin=wmin,
            wmax=wmax,
            smax=smax,
            start=start_list,
            warn_start=False,
            warn_ksub=False,
            dexp=d_exp,
            steps=num_steps,
            threads=threads,
            timeout=timeout,
            dW=dW,
            maxC=maxC,
            outC=outC_file,
            seed=seed,
            debug=debug,
            nothrottle=nothrottle,
            chunk_size=chunk_size,
            ksub=ksub,
            kwin=kwin if kwin > 0 else int(kwargs.get("win", 0) or 0),
            win_mode=win_mode,
            min_hits=min_hits,
            cov_cws=cov_cws,
            refresh=refresh,
            verbose=verbose
        )

        dist = dmin_res if (dmin_res == dmax_res or dmax_res == 0) else dmax_res
        cws = []
        if (do_cws or outC) and outC_file and os.path.exists(outC_file):
            cws = read_sparse_vectors(outC_file)
            cws.sort(key=len)
            if outC:
                _write_nzlist_file(outC, cws)

        d_info = format_bounds_list(dmin_res, dmax_res, rw_steps)

        if _use_distance_cache and code_key is not None:
            prev_steps = cached_entry.get("rw_steps", 0) if cached_entry else 0
            prev_dmax = cached_entry.get("dmax", 0) if cached_entry else 0
            prev_dmin = cached_entry.get("dmin", 0) if cached_entry else 0
            prev_cws = list(cached_entry.get("cws", [])) if cached_entry else []

            total_rw_steps = prev_steps + rw_steps
            best_dmax = (
                min(prev_dmax, dmax_res) if (prev_dmax > 0 and dmax_res > 0)
                else (dmax_res if dmax_res > 0 else prev_dmax)
            )
            best_dmin = max(prev_dmin, dmin_res)

            combined_cws = prev_cws
            if cws:
                existing_set = {tuple(cw) for cw in combined_cws}
                for cw in cws:
                    if tuple(cw) not in existing_set:
                        combined_cws.append(cw)
                        existing_set.add(tuple(cw))
                combined_cws.sort(key=len)

            d_info = format_bounds_list(best_dmin, best_dmax, total_rw_steps)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "d_info": d_info,
                "cws": combined_cws
            }
            if plain_key is not None and code_key != plain_key:
                # start-list run: merge the valid part (or all, with trust_start) into the main record
                _merge_start_into_plain(plain_key, _distance_cache[code_key], trust_start, {"": rw_steps})
            if eff_cache_file:
                save_distance_cache(eff_cache_file)

            if return_info:
                return (dist, d_info, combined_cws) if do_cws else (dist, d_info)
            return (dist, combined_cws) if do_cws else dist

        if return_info:
            return (dist, d_info, cws) if do_cws else (dist, d_info)
        return (dist, cws) if do_cws else dist

    finally:
        _remove_temp_files(temp_files, debug)


def _split_css_filename(filepath: Optional[str], sector: str) -> Optional[str]:
    """Generates sector-suffixed filename for CSS codewords (e.g. 'cws.nz' -> 'cws_X.nz')."""
    if not filepath:
        return None
    base, ext = os.path.splitext(filepath)
    if not ext:
        ext = ".nz"
    if base.endswith(f"_{sector}"):
        return f"{base}{ext}"
    elif base.endswith("_X") or base.endswith("_Z"):
        base = base[:-2]
    return f"{base}_{sector}{ext}"


def compute_css_distance(
    Hx: Any,
    Hz: Any,
    Lx: Optional[Any] = None,
    Lz: Optional[Any] = None,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmin: int = 1,
    wmax: int = 0,
    smax: Optional[int] = None,
    start: Optional[Union[int, str, Sequence[int]]] = None,
    cbeg: Optional[int] = None,
    cend: Optional[int] = None,
    noscan: int = 0,
    dW: int = -1,
    maxC: int = 0,
    finC: Optional[str] = None,
    outC: Optional[str] = None,
    do_cws: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
    solver: str = "dist_m4ri",
    codedistance_method: str = "QDistEvol",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0,
    verbose: bool = False,
    nothrottle: bool = False,
    chunk_size: int = 0,
    ksub: int = 0,
    kwin: int = 0,
    win_mode: int = 0,
    min_hits: Optional[int] = None,
    cov_cws: int = 100,
    refresh: int = 0,
    trust_start: bool = False,
    **kwargs
) -> Tuple[Any, ...]:
    """
    Computes CSS quantum code distance d = min(d_X, d_Z).

    Args:
        Hx: X-stabilizer parity check matrix.
        Hz: Z-stabilizer parity check matrix.
        Lx: Optional X-logical operator matrix (alternative to Hz as finG).
        Lz: Optional Z-logical operator matrix (alternative to Hx as finG).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance, inclusive.
        d_max / dmax: Known upper bound on distance, inclusive.
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        smax: Maximum syndrome weight for CC confinement profile.
        start: Expert option: list of CC start columns (int, comma-separated string "a,b,c", or
            sequence of ints; None or -1: all columns).  CC clusters are grown only from the listed
            columns, and they are not limited to larger column indices.  WARNING: the lower bound
            (and an exact result) is valid only if a code symmetry maps every minimum-weight codeword
            to one containing a listed column (e.g., one column per block of a quasi-cyclic code);
            the upper bound is always valid.  Results are cached under the separate key
            '<key>:start=a,b,c'; only dmax and codewords update the main cache record.
        trust_start: Expert option: accept the start-list results as valid, i.e., also merge dmin and
            rw_steps into the main cache record, and copy an existing exact start-list record to the
            main record even if no calculation is needed (default: False).
        cbeg / cend / noscan: Disabled in the Python interface: accepted for backward compatibility,
            but ignored with a warning to stderr (use the dist_m4ri binary directly).
        ksub: Experimental, should not be used (a warning is printed to stderr for ksub > 0 with method=1 or 3):
            subspace dimension sampled from ker(H) in RW (default: 0, full-matrix RW).
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        finC: Input file with initial codewords.
        outC: Output file to save codewords (NZLIST format).
        do_cws: Whether to return extracted X and Z codewords.
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug bitmap: the bits in BIN_DBG_MASK for the dist_m4ri binary, and the PY_DBG_* bits.
        verbose: Verbose reporting flag.

    Returns:
        tuple (dist, dX_info, dZ_info, cws_X, cws_Z) if do_cws
        else (dist, dX_info, dZ_info)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max
    start_list = _prepare_start_option(start, trust_start, noscan, cbeg, cend, solver)
    _warn_experimental_ksub(ksub, method, solver)

    can_compute_Z = (
        Hx is not None
        and (hasattr(Hx, 'shape') and Hx.shape[0] > 0 if not isinstance(Hx, (str, Path)) else True)
    )
    can_compute_X = (
        Hz is not None
        and (hasattr(Hz, 'shape') and Hz.shape[0] > 0 if not isinstance(Hz, (str, Path)) else True)
    )

    if not can_compute_Z and not can_compute_X:
        raise ValueError("Cannot compute CSS distance: Both Hx and Hz are empty.")
    css_sectors = [s for s, ok in (("X", can_compute_X), ("Z", can_compute_Z)) if ok]

    finC = check_finc_outc(finC, outC, verbose=verbose)

    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None
    plain_key = None
    plain_entry = None

    if solver == "codedistance":
        if do_cws or outC:
            raise ValueError("Codeword extraction is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
        import numpy as np
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        dist_Z, dist_X = None, None
        dX_info, dZ_info = None, None

        Hx_mat = Hx.toarray() if hasattr(Hx, 'toarray') else (
            np.asarray(Hx, dtype=np.int8) if isinstance(Hx, (np.ndarray, list)) else None
        )
        Hz_mat = Hz.toarray() if hasattr(Hz, 'toarray') else (
            np.asarray(Hz, dtype=np.int8) if isinstance(Hz, (np.ndarray, list)) else None
        )

        if can_compute_Z:
            res_Z = codedistance.CSScodeDistance(
                Hx_mat, Hz_mat, method=codedistance_method, params=params.copy(),
                component="Z", seed=seed if seed != 0 else None
            )
            dist_Z = res_Z.get("d", -1)
            dZ_info = format_bounds_list(dist_Z, dist_Z, 0) if dist_Z > 0 else [0, 0, 0]

        if can_compute_X:
            res_X = codedistance.CSScodeDistance(
                Hx_mat, Hz_mat, method=codedistance_method, params=params.copy(),
                component="X", seed=seed if seed != 0 else None
            )
            dist_X = res_X.get("d", -1)
            dX_info = format_bounds_list(dist_X, dist_X, 0) if dist_X > 0 else [0, 0, 0]

        if dist_X is not None and dist_Z is not None:
            dist = min(dist_Z, dist_X) if (dist_Z > 0 and dist_X > 0) else max(dist_Z, dist_X)
        elif dist_X is not None:
            dist = dist_X
        else:
            dist = dist_Z

        return (dist, dX_info, dZ_info)

    # Solver is native multithreaded dist_m4ri
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            hx_state = get_sparse_array_state(Hx) if can_compute_Z else "none"
            hz_state = get_sparse_array_state(Hz) if can_compute_X else "none"
            plain_key = f"css:X={hx_state}:Z={hz_state}"
            if Lx is not None or Lz is not None:
                lx_state = get_sparse_array_state(Lx) if Lx is not None else "none"
                lz_state = get_sparse_array_state(Lz) if Lz is not None else "none"
                plain_key = f"{plain_key}:Lx={lx_state}:Lz={lz_state}"
            code_key, cached_entry, plain_entry = _resolve_start_cache(
                plain_key, start_list,
                lambda e: _css_entry_is_exact(e, can_compute_X, can_compute_Z, bool(do_cws or outC)), trust_start,
                lambda e: _merge_start_into_plain(plain_key, e, True, None, css_sectors), eff_cache_file, verbose
            )
            if cached_entry is not None:
                c_dmin_x = cached_entry.get("dmin_X", cached_entry.get("dmin", 0))
                c_dmax_x = cached_entry.get("dmax_X", cached_entry.get("dmax", 0))
                c_dmin_z = cached_entry.get("dmin_Z", cached_entry.get("dmin", 0))
                c_dmax_z = cached_entry.get("dmax_Z", cached_entry.get("dmax", 0))
                x_exact = (not can_compute_X) or (c_dmin_x > 0 and c_dmin_x == c_dmax_x)
                z_exact = (not can_compute_Z) or (c_dmin_z > 0 and c_dmin_z == c_dmax_z)
                # If exact distance is already proven for all requested sectors and not asking for more codewords
                if x_exact and z_exact and (c_dmax_x > 0 or c_dmax_z > 0):
                    if not (do_cws or outC) or (cached_entry.get("cws_X") and cached_entry.get("cws_Z")):
                        dx_res = cached_entry.get(
                            "dX",
                            format_bounds_list(
                                c_dmin_x,
                                c_dmax_x,
                                cached_entry.get("rw_steps_X", 0)
                            )
                        )
                        dz_res = cached_entry.get(
                            "dZ",
                            format_bounds_list(
                                c_dmin_z,
                                c_dmax_z,
                                cached_entry.get("rw_steps_Z", 0)
                            )
                        )
                        if verbose:
                            print(f"[dist_m4ri] Cache retrieval: SUCCESS "
                                  f"(found cached exact CSS distance for '{code_key}')")
                            print(f"[dist_m4ri] Cached result: dist={cached_entry['dist']}, "
                                  f"dX={format_bounds_str(dx_res)}, dZ={format_bounds_str(dz_res)}")
                        elif debug & PY_DBG_CACHE:
                            print("[dist_m4ri] Cache hit for CSS distance (exact distance known)!")
                        cws_x = cached_entry.get("cws_X", [])
                        cws_z = cached_entry.get("cws_Z", [])
                        if outC:
                            existing_cws = (
                                read_sparse_vectors(finC)
                                if (finC and os.path.exists(finC)
                                    and (finC == outC or os.path.abspath(finC) == os.path.abspath(outC)))
                                else []
                            )
                            combined = existing_cws + (cws_x or []) + (cws_z or [])
                            seen = set()
                            unique = []
                            for cw in combined:
                                t = tuple(cw)
                                if t not in seen:
                                    seen.add(t)
                                    unique.append(cw)
                            unique.sort(key=len)
                            if unique:
                                _write_nzlist_file(outC, unique)
                        return (
                            cached_entry["dist"], dx_res, dz_res,
                            cws_x, cws_z
                        ) if do_cws else (
                            cached_entry["dist"], dx_res, dz_res
                        )
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: PARTIAL (cached CSS bounds: "
                          f"dmin={cached_entry.get('dmin', 0)}, dmax={cached_entry.get('dmax', 0)}; "
                          f"continuing search)")
            else:
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: MISS (no entry for '{code_key}')")
        except Exception:
            code_key = None
            cached_entry = None
            plain_key = None
            plain_entry = None
    else:
        if verbose:
            print("[dist_m4ri] Cache retrieval: DISABLED (cache is turned off)")

    temp_files = []
    try:
        if can_compute_Z and not isinstance(Hx, (str, Path)):
            file_Hx = _matrix_to_file(Hx, extension="_Hx.mtx")
            temp_files.append(file_Hx)
        else:
            file_Hx = str(Hx) if can_compute_Z else None

        if can_compute_X and not isinstance(Hz, (str, Path)):
            file_Hz = _matrix_to_file(Hz, extension="_Hz.mtx")
            temp_files.append(file_Hz)
        else:
            file_Hz = str(Hz) if can_compute_X else None

        if Lx is not None and not isinstance(Lx, (str, Path)):
            file_Lx = _matrix_to_file(Lx, extension="_Lx.mtx")
            temp_files.append(file_Lx)
        else:
            file_Lx = str(Lx) if Lx is not None else None

        if Lz is not None and not isinstance(Lz, (str, Path)):
            file_Lz = _matrix_to_file(Lz, extension="_Lz.mtx")
            temp_files.append(file_Lz)
        else:
            file_Lz = str(Lz) if Lz is not None else None

        outZ = create_unique_file(extension="_Z.nz") if ((do_cws or outC) and can_compute_Z) else None
        outX = create_unique_file(extension="_X.nz") if ((do_cws or outC) and can_compute_X) else None
        if outZ:
            temp_files.append(outZ)
        if outX:
            temp_files.append(outX)

        dist_Z, dist_X = None, None
        dmin_z, dmax_z, rw_steps_z = 0, 0, 0
        dmin_x, dmax_x, rw_steps_x = 0, 0, 0
        cws_Z, cws_X = [], []

        # Resolve sector-specific input codewords (finC_Z and finC_X)
        finC_Z = None
        finC_X = None
        if finC:
            outC_Z_name = _split_css_filename(outC, "Z") if outC else None
            outC_X_name = _split_css_filename(outC, "X") if outC else None
            cand_Z = _split_css_filename(finC, "Z")
            cand_X = _split_css_filename(finC, "X")

            finC_Z = check_finc_outc(cand_Z, outC_Z_name, verbose=False)
            finC_X = check_finc_outc(cand_X, outC_X_name, verbose=False)

            # Fallback: if sector-suffixed files don't exist, check raw finC
            if not finC_Z and not finC_X and os.path.exists(finC):
                finC_Z = check_finc_outc(finC, outC_Z_name, verbose=False)
                finC_X = check_finc_outc(finC, outC_X_name, verbose=False)

        # Seed sector-specific bounds from cache if available
        eff_dmin_z, eff_dmax_z = eff_dmin, eff_dmax
        eff_dmin_x, eff_dmax_x = eff_dmin, eff_dmax
        if cached_entry is not None:
            cz_min = cached_entry.get("dmin_Z", cached_entry.get("dmin", 0))
            cz_max = cached_entry.get("dmax_Z", cached_entry.get("dmax", 0))
            cx_min = cached_entry.get("dmin_X", cached_entry.get("dmin", 0))
            cx_max = cached_entry.get("dmax_X", cached_entry.get("dmax", 0))
            if cz_min > 1:
                eff_dmin_z = max(eff_dmin_z, cz_min) if eff_dmin_z > 1 else cz_min
            if cz_max > 0:
                eff_dmax_z = min(eff_dmax_z, cz_max) if eff_dmax_z > 0 else cz_max
            if cx_min > 1:
                eff_dmin_x = max(eff_dmin_x, cx_min) if eff_dmin_x > 1 else cx_min
            if cx_max > 0:
                eff_dmax_x = min(eff_dmax_x, cx_max) if eff_dmax_x > 0 else cx_max
        # Start-list run: also seed from the main (plain-key) record
        eff_dmin_z, eff_dmax_z = _seed_bounds_from_entry(eff_dmin_z, eff_dmax_z, plain_entry, "_Z")
        eff_dmin_x, eff_dmax_x = _seed_bounds_from_entry(eff_dmin_x, eff_dmax_x, plain_entry, "_X")

        # Z-distance: Hx as finH, Hz as finG (or Lx as finL dual logical operators)
        if can_compute_Z:
            dmin_z, dmax_z, rw_steps_z = run_dist_m4ri(
                dist_m4ri_path=dist_m4ri,
                method=method,
                finH=file_Hx,
                finG=file_Hz if file_Lx is None else None,
                finL=file_Lx,
                finC=finC_Z,
                dmin=eff_dmin_z,
                dmax=eff_dmax_z,
                wmin=wmin,
                wmax=wmax,
                smax=smax,
                start=start_list,
                warn_start=False,
                warn_ksub=False,
                dexp=d_exp,
                steps=num_steps,
                threads=threads,
                timeout=timeout,
                dW=dW,
                maxC=maxC,
                outC=outZ,
                seed=seed,
                debug=debug,
                nothrottle=nothrottle,
                chunk_size=chunk_size,
                ksub=ksub,
                kwin=kwin if kwin > 0 else int(kwargs.get("win", 0) or 0),
                win_mode=win_mode,
                min_hits=min_hits,
                cov_cws=cov_cws,
                refresh=refresh,
                verbose=verbose
            )
            _last_css_stats["Z"] = dict(_last_run_stats)
            dist_Z = dmin_z if (dmin_z == dmax_z or dmax_z == 0) else dmax_z
            if (do_cws or outC) and outZ and os.path.exists(outZ):
                cws_Z = read_sparse_vectors(outZ)
                cws_Z.sort(key=len)

        # X-distance: Hz as finH, Hx as finG (or Lz as finL dual logical operators)
        if can_compute_X:
            dmin_x, dmax_x, rw_steps_x = run_dist_m4ri(
                dist_m4ri_path=dist_m4ri,
                method=method,
                finH=file_Hz,
                finG=file_Hx if file_Lz is None else None,
                finL=file_Lz,
                finC=finC_X,
                dmin=eff_dmin_x,
                dmax=eff_dmax_x,
                wmin=wmin,
                wmax=wmax,
                smax=smax,
                start=start_list,
                warn_start=False,
                warn_ksub=False,
                dexp=d_exp,
                steps=num_steps,
                threads=threads,
                timeout=timeout,
                dW=dW,
                maxC=maxC,
                outC=outX,
                seed=seed,
                debug=debug,
                nothrottle=nothrottle,
                chunk_size=chunk_size,
                ksub=ksub,
                kwin=kwin if kwin > 0 else int(kwargs.get("win", 0) or 0),
                win_mode=win_mode,
                min_hits=min_hits,
                cov_cws=cov_cws,
                refresh=refresh,
                verbose=verbose
            )
            _last_css_stats["X"] = dict(_last_run_stats)
            dist_X = dmin_x if (dmin_x == dmax_x or dmax_x == 0) else dmax_x
            if (do_cws or outC) and outX and os.path.exists(outX):
                cws_X = read_sparse_vectors(outX)
                cws_X.sort(key=len)

        run_rw_steps = {"_X": rw_steps_x, "_Z": rw_steps_z}  # RW steps of this run (before merging the cache)
        if _use_distance_cache and cached_entry is not None:
            if can_compute_Z:
                pz_min = cached_entry.get("dmin_Z", 0)
                pz_max = cached_entry.get("dmax_Z", 0)
                pz_rw = cached_entry.get("rw_steps_Z", 0)
                dmin_z = max(pz_min, dmin_z)
                dmax_z = min(pz_max, dmax_z) if (pz_max > 0 and dmax_z > 0) else (dmax_z if dmax_z > 0 else pz_max)
                if dmax_z > 0 and dmin_z >= dmax_z:
                    dmin_z = dmax_z
                rw_steps_z = 0 if (dmin_z > 0 and dmin_z == dmax_z) else (pz_rw + rw_steps_z)
                dist_Z = dmin_z if (dmin_z == dmax_z or dmax_z == 0) else dmax_z
            if can_compute_X:
                px_min = cached_entry.get("dmin_X", 0)
                px_max = cached_entry.get("dmax_X", 0)
                px_rw = cached_entry.get("rw_steps_X", 0)
                dmin_x = max(px_min, dmin_x)
                dmax_x = min(px_max, dmax_x) if (px_max > 0 and dmax_x > 0) else (dmax_x if dmax_x > 0 else px_max)
                if dmax_x > 0 and dmin_x >= dmax_x:
                    dmin_x = dmax_x
                rw_steps_x = 0 if (dmin_x > 0 and dmin_x == dmax_x) else (px_rw + rw_steps_x)
                dist_X = dmin_x if (dmin_x == dmax_x or dmax_x == 0) else dmax_x

        dX_info = format_bounds_list(dmin_x, dmax_x, rw_steps_x) if can_compute_X else None
        dZ_info = format_bounds_list(dmin_z, dmax_z, rw_steps_z) if can_compute_Z else None

        if dist_X is not None and dist_Z is not None:
            dist = min(dist_Z, dist_X) if (dist_Z > 0 and dist_X > 0) else max(dist_Z, dist_X)
        elif dist_X is not None:
            dist = dist_X
        else:
            dist = dist_Z

        if outC:
            outC_X = _split_css_filename(outC, "X")
            outC_Z = _split_css_filename(outC, "Z")
            if cws_X:
                _write_nzlist_file(outC_X, cws_X)
            if cws_Z:
                _write_nzlist_file(outC_Z, cws_Z)

        res_tuple = (dist, dX_info, dZ_info, cws_X, cws_Z) if do_cws else (dist, dX_info, dZ_info)

        if _use_distance_cache and code_key is not None:
            total_rw_steps = (rw_steps_z if can_compute_Z else 0) + (rw_steps_x if can_compute_X else 0)

            curr_dmax = 0
            if can_compute_Z and can_compute_X:
                if dmax_z > 0 and dmax_x > 0:
                    curr_dmax = min(dmax_z, dmax_x)
            elif can_compute_Z:
                curr_dmax = dmax_z
            elif can_compute_X:
                curr_dmax = dmax_x
            best_dmax = curr_dmax

            curr_dmin = 0
            if can_compute_Z and can_compute_X:
                curr_dmin = min(dmin_z, dmin_x)
            elif can_compute_Z:
                curr_dmin = dmin_z
            elif can_compute_X:
                curr_dmin = dmin_x
            best_dmin = curr_dmin

            combined_cws_x = list(cached_entry.get("cws_X", [])) if cached_entry else []
            if cws_X:
                existing_x = {tuple(cw) for cw in combined_cws_x}
                for cw in cws_X:
                    if tuple(cw) not in existing_x:
                        combined_cws_x.append(cw)
                        existing_x.add(tuple(cw))
                combined_cws_x.sort(key=len)

            combined_cws_z = list(cached_entry.get("cws_Z", [])) if cached_entry else []
            if cws_Z:
                existing_z = {tuple(cw) for cw in combined_cws_z}
                for cw in cws_Z:
                    if tuple(cw) not in existing_z:
                        combined_cws_z.append(cw)
                        existing_z.add(tuple(cw))
                combined_cws_z.sort(key=len)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "dmin_X": dmin_x,
                "dmax_X": dmax_x,
                "rw_steps_X": rw_steps_x,
                "dmin_Z": dmin_z,
                "dmax_Z": dmax_z,
                "rw_steps_Z": rw_steps_z,
                "dX": dX_info,
                "dZ": dZ_info,
                "cws_X": combined_cws_x,
                "cws_Z": combined_cws_z
            }
            if plain_key is not None and code_key != plain_key:
                # start-list run: merge the valid part (or all, with trust_start) into the main record
                _merge_start_into_plain(plain_key, _distance_cache[code_key], trust_start, run_rw_steps, css_sectors)
            if eff_cache_file:
                save_distance_cache(eff_cache_file)
        return res_tuple

    finally:
        _remove_temp_files(temp_files, debug)


def _resolve_stim_dem_out_path(
    out_spec: Optional[Union[bool, str, Path]],
    circuit_src: str,
    mode_tag: str,
    ext: str,
    out_dir: Optional[Union[str, Path]] = None
) -> Optional[str]:
    """Resolves output filepath for --out-dem or --out-stim."""
    if out_spec is None or out_spec is False:
        return None
    if isinstance(out_spec, str) and out_spec.strip().lower() in (
        "0", "false", "no", "off", "none"
    ):
        return None

    src_path = (
        Path(circuit_src)
        if (circuit_src and circuit_src != "stim.Circuit")
        else None
    )
    eff_out_dir = (
        Path(out_dir)
        if out_dir is not None
        else (src_path.parent if src_path is not None else Path("."))
    )

    is_auto = (out_spec is True) or (
        isinstance(out_spec, str)
        and out_spec.strip().lower() in ("", "1", "true", "yes", "auto")
    )
    if is_auto:
        if src_path is not None:
            stem = src_path.stem
            if stem.endswith("_full") or stem.endswith("_simp"):
                stem = stem[:-5]
        else:
            stem = "circuit"
        out_path = eff_out_dir / f"{stem}{mode_tag}{ext}"
    else:
        spec_str = str(out_spec)
        spec_path = Path(spec_str)
        if not spec_path.is_absolute() and os.path.dirname(spec_str) == "":
            out_path = eff_out_dir / spec_path
        else:
            out_path = spec_path

    parent_dir = os.path.dirname(os.path.abspath(str(out_path)))
    if parent_dir:
        os.makedirs(parent_dir, exist_ok=True)
    return str(out_path)


def compute_dem_distance(
    dem: Optional[Any] = None,
    circuit: Optional[Any] = None,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmin: int = 1,
    wmax: int = 0,
    smax: Optional[int] = None,
    start: Optional[Union[int, str, Sequence[int]]] = None,
    cbeg: Optional[int] = None,
    cend: Optional[int] = None,
    noscan: int = 0,
    dW: int = -1,
    maxC: int = 0,
    pmin: float = 0.0,
    finC: Optional[str] = None,
    outC: Optional[str] = None,
    do_cws: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
    solver: str = "dist_m4ri",
    codedistance_method: str = "UndetectableErrorStim",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0,
    verbose: bool = False,
    nothrottle: bool = False,
    chunk_size: int = 0,
    ksub: int = 0,
    kwin: int = 0,
    win_mode: int = 0,
    min_hits: Optional[int] = None,
    cov_cws: int = 100,
    refresh: int = 0,
    simple: Optional[bool] = None,
    full: bool = False,
    basis: Optional[str] = None,
    rounds: Optional[int] = None,
    out_dir: Optional[Union[str, Path]] = None,
    out_dem: Optional[Union[bool, str, Path]] = None,
    out_stim: Optional[Union[bool, str, Path]] = None,
    trust_start: bool = False,
    **kwargs
) -> Tuple[Any, ...]:
    """
    Computes minimum graph/hypergraph distance of a Stim Detector Error Model (DEM).

    Args:
        dem: stim.DetectorErrorModel object or path to .dem / .stim file.
        circuit: stim.Circuit object or path to .stim file (converted to DEM).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance, inclusive.
        d_max / dmax: Known upper bound on distance, inclusive.
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        smax: Maximum syndrome weight for CC confinement profile.
        start: Expert option: list of CC start columns (int, comma-separated string "a,b,c", or
            sequence of ints; None or -1: all columns).  CC clusters are grown only from the listed
            columns, and they are not limited to larger column indices.  WARNING: the lower bound
            (and an exact result) is valid only if a code symmetry maps every minimum-weight codeword
            to one containing a listed column (e.g., one column per block of a quasi-cyclic code);
            the upper bound is always valid.  Results are cached under the separate key
            '<key>:start=a,b,c'; only dmax and codewords update the main cache record.
        trust_start: Expert option: accept the start-list results as valid, i.e., also merge dmin and
            rw_steps into the main cache record, and copy an existing exact start-list record to the
            main record even if no calculation is needed (default: False).
        cbeg / cend / noscan: Disabled in the Python interface: accepted for backward compatibility,
            but ignored with a warning to stderr (use the dist_m4ri binary directly).
        ksub: Experimental, should not be used (a warning is printed to stderr for ksub > 0 with method=1 or 3):
            subspace dimension sampled from ker(H) in RW (default: 0, full-matrix RW).
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        pmin: Probability cutoff for error mechanisms in DEM.
            Note: CC start columns (start) index the error mechanisms in their order of appearance in
            the (flattened) DEM, after the pmin cutoff.
        finC: Input file with initial codewords.
        outC: Output file to save codewords (NZLIST format).
        do_cws: Whether to return extracted error mechanisms / codewords.
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug bitmap: the bits in BIN_DBG_MASK for the dist_m4ri binary, and the PY_DBG_* bits.
        verbose: Verbose reporting flag.
        simple: If True (default for CSS .stim circuits), keep only primary-basis detectors.
        full: If True, keep all detectors in .stim circuits.
        basis: Optional primary memory basis override ('X' or 'Z') for .stim circuits.
        rounds: Optional number of REPEAT rounds for .stim circuits (warns if no REPEAT block).
        out_dir: Optional output directory for saved DEM and/or Stim circuit files.
        out_dem: Optional bool or filename to save the constructed DEM (True = auto-named
            <basename>_simp.dem or <basename>_full.dem).
        out_stim: Optional bool or filename to save the processed (noisy) Stim circuit
            (True = auto-named <basename>_simp.stim or <basename>_full.stim).

    Returns:
        tuple (dist, d_info, cws) if do_cws else (dist, d_info)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max
    start_list = _prepare_start_option(start, trust_start, noscan, cbeg, cend, solver)
    _warn_experimental_ksub(ksub, method, solver)

    if dem is not None and isinstance(dem, (str, Path)) and str(dem).endswith(".stim"):
        circuit = dem
        dem = None

    if rounds is not None and circuit is None:
        raise ValueError(
            "--rounds option is only supported for Stim circuit (.stim) inputs."
        )

    saved_dem_path: Optional[str] = None
    circ_verbose = verbose or bool(debug & PY_DBG_CIRCUITS)  # Stim circuit to DEM conversion details

    if dem is None and circuit is not None:
        circuit_src = str(circuit) if isinstance(circuit, (str, Path)) else "stim.Circuit"
        if isinstance(circuit, (str, Path)):
            stim = _get_stim()
            circuit = stim.Circuit.from_file(str(circuit))
        if hasattr(circuit, 'detector_error_model'):
            has_rep = False
            if rounds is not None:
                if rounds < 0:
                    raise ValueError(f"Invalid rounds={rounds}; must be >= 0.")
                has_rep = has_repeat_block(circuit)
                if not has_rep:
                    sys.stderr.write(
                        f"# Warning: Circuit '{circuit_src}' does not contain a "
                        f"REPEAT block; ignoring rounds={rounds}.\n"
                    )

            circuit, empty_removed = remove_empty_detectors(circuit)
            if empty_removed > 0 and circ_verbose:
                print(
                    f"[dist_m4ri] Removed {empty_removed} empty DETECTOR(s) "
                    f"from '{circuit_src}'"
                )

            filepath_hint = circuit_src if circuit_src != "stim.Circuit" else None
            t_res = classify_qubits_thorough(
                circuit,
                verbose=False,
                basis_arg=basis,
                filepath=filepath_hint
            )
            eff_basis = t_res["basis"]
            is_css_circ = t_res["is_css"]
            is_rot_css = t_res.get("is_rotated_css", False)

            if rounds is not None and has_rep:
                circuit = set_circuit_rounds(circuit, rounds)

            if full:
                use_simple = False
            elif simple is not None:
                use_simple = bool(simple)
            else:
                use_simple = is_css_circ

            stripped_det = 0
            kept_det = circuit.num_detectors
            if use_simple and is_css_circ:
                circuit, stripped_det, kept_det = strip_minority_detectors(
                    circuit,
                    eff_basis,
                    verbose=circ_verbose,
                    thorough=True,
                    thorough_res=t_res
                )
            elif simple is True and not is_css_circ and circ_verbose:
                print(
                    f"[dist_m4ri] Warning: --simple requested, but '{circuit_src}' "
                    "is not classified as a CSS circuit; keeping all detectors."
                )

            mode_tag = "_simp" if (use_simple and is_css_circ) else "_full"

            if circ_verbose:
                mode_str = (
                    "simple (primary-basis detectors only)"
                    if (use_simple and is_css_circ)
                    else "full (all detectors)"
                )
                css_type = (
                    "rotated-CSS (e.g. XZZX)"
                    if is_rot_css
                    else ("CSS" if is_css_circ else "non-CSS")
                )
                rot_info = (
                    f", rotated_data={len(t_res.get('rotated_data_qubits', []))}"
                    if is_rot_css
                    else ""
                )
                rnd_cnt = count_circuit_rounds(circuit)
                print(
                    f"[dist_m4ri] Circuit basis tracking: type={css_type}, "
                    f"basis={eff_basis}, data={len(t_res['data_qubits'])}, "
                    f"X_anc={len(t_res['x_ancillas'])}, "
                    f"Z_anc={len(t_res['z_ancillas'])}, "
                    f"routing={len(t_res['routing_qubits'])}{rot_info}, "
                    f"rounds={rnd_cnt}"
                )
                print(
                    f"[dist_m4ri] Circuit mode: {mode_str} "
                    f"(kept={kept_det}, stripped={stripped_det} detectors)"
                )

            noise_added = False
            p_noise = float(kwargs.get("p_noise", 0.001))
            if not has_noise(circuit):
                circuit = add_noise(circuit, p=p_noise)
                noise_added = True

            stim_out_path = _resolve_stim_dem_out_path(
                out_stim, circuit_src, mode_tag, ".stim", out_dir=out_dir
            )
            if stim_out_path is not None:
                circuit.to_file(stim_out_path)
                if circ_verbose:
                    print(f"[dist_m4ri] Saved Stim circuit to '{stim_out_path}'")

            try:
                dem = circuit.detector_error_model(decompose_errors=True)
                decomp_used = True
            except Exception:
                dem = circuit.detector_error_model(decompose_errors=False)
                decomp_used = False
            if circ_verbose:
                noise_msg = (
                    f"added circuit-level noise (p={p_noise})"
                    if noise_added else "existing noise detected"
                )
                print(
                    f"[dist_m4ri] Converted '{circuit_src}' to DEM "
                    f"({noise_msg}, decompose_errors={decomp_used})"
                )
                if all(hasattr(dem, a) for a in ("num_detectors", "num_observables", "num_errors")):
                    print(
                        f"[dist_m4ri] DEM size: {dem.num_detectors} detectors, {dem.num_observables} observables, "
                        f"{dem.num_errors} error mechanisms"
                    )

            dem_out_path = _resolve_stim_dem_out_path(
                out_dem, circuit_src, mode_tag, ".dem", out_dir=out_dir
            )
            if dem_out_path is not None:
                if hasattr(dem, 'flattened'):
                    dem.flattened().to_file(dem_out_path)
                elif hasattr(dem, 'to_file'):
                    dem.to_file(dem_out_path)
                else:
                    with open(dem_out_path, 'w') as f_dem:
                        f_dem.write(str(dem))
                saved_dem_path = dem_out_path
                if circ_verbose:
                    print(f"[dist_m4ri] Saved DEM to '{dem_out_path}'")
        else:
            raise ValueError("Provided circuit object does not have detector_error_model() method.")
    elif dem is not None and out_dem is not None and out_dem is not False:
        dem_src = str(dem) if isinstance(dem, (str, Path)) else "stim.Circuit"
        dem_out_path = _resolve_stim_dem_out_path(
            out_dem, dem_src, "_full", ".dem", out_dir=out_dir
        )
        if dem_out_path is not None:
            if isinstance(dem, (str, Path)) and os.path.exists(str(dem)):
                if os.path.abspath(str(dem)) != os.path.abspath(dem_out_path):
                    shutil.copyfile(str(dem), dem_out_path)
            elif hasattr(dem, 'flattened'):
                dem.flattened().to_file(dem_out_path)
            elif hasattr(dem, 'to_file'):
                dem.to_file(dem_out_path)
            else:
                with open(dem_out_path, 'w') as f_dem:
                    f_dem.write(str(dem))
            saved_dem_path = dem_out_path
            if circ_verbose:
                print(f"[dist_m4ri] Saved DEM to '{dem_out_path}'")

    if dem is None:
        raise ValueError("Either 'dem' or 'circuit' must be provided.")

    finC = check_finc_outc(finC, outC, verbose=verbose)

    if solver == "codedistance":
        if do_cws or outC:
            raise ValueError("Codeword extraction is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
        params = dict(codedistance_params or {})
        params.setdefault("filterCircuit", False)
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        if circuit is not None:
            res = codedistance.circuitDistance(
                circuit, method=codedistance_method, params=params,
                seed=seed if seed != 0 else None
            )
        else:
            H, L, priors = codedistance.StimDEM2HL(dem)
            if "priors" not in params and len(priors) > 0:
                params["priors"] = priors
            res = codedistance.codeDistance(
                H, L, tB=1, method=codedistance_method, params=params,
                seed=seed if seed != 0 else None
            )
        d = res.get("d", -1)
        d_info = format_bounds_list(d, d, 0) if d > 0 else [0, 0, 0]
        return d, d_info

    # Solver is native multithreaded dist_m4ri
    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None
    plain_key = None
    plain_entry = None
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            dem_obj = dem if dem is not None else circuit
            dem_state = get_sparse_array_state(dem_obj)
            plain_key = f"dem:{dem_state}" if pmin <= 0.0 else f"dem:{dem_state}:pmin={pmin}"
            code_key, cached_entry, plain_entry = _resolve_start_cache(
                plain_key, start_list, lambda e: _single_entry_is_exact(e, bool(do_cws or outC)), trust_start,
                lambda e: _merge_start_into_plain(plain_key, e, True), eff_cache_file, verbose
            )
            eff_dmin, eff_dmax = _seed_bounds_from_entry(eff_dmin, eff_dmax, plain_entry)
            if cached_entry is not None:
                # If exact distance is already proven and not asking for more codewords
                if cached_entry.get("dmin", 0) > 0 and cached_entry.get("dmin") == cached_entry.get("dmax"):
                    if not (do_cws or outC) or (cached_entry.get("cws") and len(cached_entry["cws"]) > 0):
                        d_info = format_bounds_list(
                            cached_entry.get("dmin", 0), cached_entry.get("dmax", 0), cached_entry.get("rw_steps", 0)
                        )
                        if verbose:
                            print(f"[dist_m4ri] Cache retrieval: SUCCESS "
                                  f"(found cached exact DEM distance for '{code_key}')")
                            print(f"[dist_m4ri] Cached result: dist={cached_entry['dist']}, "
                                  f"bounds={format_bounds_str(d_info)}")
                        elif debug & PY_DBG_CACHE:
                            print("[dist_m4ri] Cache hit for DEM distance (exact distance known)!")
                        cws_res = cached_entry.get("cws", [])
                        if outC and cws_res:
                            _write_nzlist_file(outC, cws_res)
                        return (
                            cached_entry["dist"], d_info, cws_res
                        ) if do_cws else (
                            cached_entry["dist"], d_info
                        )
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: PARTIAL (cached DEM bounds: "
                          f"dmin={cached_entry.get('dmin', 0)}, dmax={cached_entry.get('dmax', 0)}; "
                          f"continuing search)")
                # Seed bounds from cache
                if eff_dmax == 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = cached_entry["dmax"]
                elif eff_dmax > 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = min(eff_dmax, cached_entry["dmax"])
                if eff_dmin <= 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = cached_entry["dmin"]
                elif eff_dmin > 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = max(eff_dmin, cached_entry["dmin"])
            else:
                if verbose:
                    print(f"[dist_m4ri] Cache retrieval: MISS (no entry for '{code_key}')")
        except Exception:
            code_key = None
            cached_entry = None
            plain_key = None
            plain_entry = None
    else:
        if verbose:
            print("[dist_m4ri] Cache retrieval: DISABLED (cache is turned off)")

    temp_files = []
    try:
        if saved_dem_path is not None and os.path.exists(saved_dem_path):
            file_dem = saved_dem_path
        elif isinstance(dem, (str, Path)) and os.path.exists(str(dem)):
            file_dem = str(dem)
        else:
            file_dem = create_unique_file(extension=".dem")
            temp_files.append(file_dem)
            if hasattr(dem, 'flattened'):
                dem.flattened().to_file(file_dem)
            elif hasattr(dem, 'to_file'):
                dem.to_file(file_dem)
            else:
                with open(file_dem, 'w') as f:
                    f.write(str(dem))

        outC_file = None
        if do_cws or outC:
            outC_file = create_unique_file(extension="_out.nz")
            temp_files.append(outC_file)

        dmin_res, dmax_res, rw_steps = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            fdem=file_dem,
            finC=finC,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmin=wmin,
            wmax=wmax,
            smax=smax if smax is not None else 0,
            start=start_list,
            warn_start=False,
            warn_ksub=False,
            dexp=d_exp,
            steps=num_steps,
            threads=threads,
            timeout=timeout,
            pmin=pmin,
            dW=dW,
            maxC=maxC,
            outC=outC_file,
            seed=seed,
            debug=debug,
            nothrottle=nothrottle,
            chunk_size=chunk_size,
            ksub=ksub,
            kwin=kwin if kwin > 0 else int(kwargs.get("win", 0) or 0),
            win_mode=win_mode,
            min_hits=min_hits,
            cov_cws=cov_cws,
            refresh=refresh,
            verbose=verbose
        )

        dist = dmin_res if (dmin_res == dmax_res or dmax_res == 0) else dmax_res
        cws = []
        if (do_cws or outC) and outC_file and os.path.exists(outC_file):
            cws = read_sparse_vectors(outC_file)
            cws.sort(key=len)
            if outC:
                _write_nzlist_file(outC, cws)

        d_info = format_bounds_list(dmin_res, dmax_res, rw_steps)

        if _use_distance_cache and code_key is not None:
            prev_steps = cached_entry.get("rw_steps", 0) if cached_entry else 0
            prev_dmax = cached_entry.get("dmax", 0) if cached_entry else 0
            prev_dmin = cached_entry.get("dmin", 0) if cached_entry else 0
            prev_cws = list(cached_entry.get("cws", [])) if cached_entry else []

            total_rw_steps = prev_steps + rw_steps
            best_dmax = (
                min(prev_dmax, dmax_res) if (prev_dmax > 0 and dmax_res > 0)
                else (dmax_res if dmax_res > 0 else prev_dmax)
            )
            best_dmin = max(prev_dmin, dmin_res)

            combined_cws = prev_cws
            if cws:
                existing_set = {tuple(cw) for cw in combined_cws}
                for cw in cws:
                    if tuple(cw) not in existing_set:
                        combined_cws.append(cw)
                        existing_set.add(tuple(cw))
                combined_cws.sort(key=len)

            d_info = format_bounds_list(best_dmin, best_dmax, total_rw_steps)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "d_info": d_info,
                "cws": combined_cws
            }
            if plain_key is not None and code_key != plain_key:
                # start-list run: merge the valid part (or all, with trust_start) into the main record
                _merge_start_into_plain(plain_key, _distance_cache[code_key], trust_start, {"": rw_steps})
            if eff_cache_file:
                save_distance_cache(eff_cache_file)

        if do_cws:
            return dist, d_info, cws
        return dist, d_info

    finally:
        _remove_temp_files(temp_files, debug)


def _write_nzlist_file(filepath: str, cws: List[List[int]]) -> None:
    """Writes codewords to a text file in NZLIST format (1-based indices).
    Skips creating the file if cws is empty.
    """
    if not cws:
        return
    with open(filepath, "w") as f:
        f.write("%% NZLIST\n")
        f.write(f"% {len(cws)} codewords\n")
        for cw in cws:
            f.write(f"{len(cw)} " + " ".join(str(idx + 1) for idx in cw) + "\n")


def parse_cli_args(argv: List[str]) -> Dict[str, Any]:
    """Parses CLI arguments supporting both key=value pairs and standard --flags."""
    args: Dict[str, Any] = {
        "method": 3,
        "finH": None,
        "finG": None,
        "finL": None,
        "fin": None,
        "fdem": None,
        "Hx": None,
        "Hz": None,
        "Lx": None,
        "Lz": None,
        "finC": None,
        "outC": None,
        "dmin": 0,
        "dmax": 0,
        "wmin": 1,
        "wmax": 0,
        "smax": None,
        "start": None,
        "cbeg": None,
        "cend": None,
        "trust_start": False,
        "css": None,
        "dexp": 0,
        "steps": None,
        "threads": None,
        "timeout": 60.0,
        "dW": -1,
        "maxC": 0,
        "pmin": 0.0,
        "noscan": 0,
        "classical": -1,
        "seed": 0,
        "debug": 0,
        "solver": "dist_m4ri",
        "cache_file": "tmp_dist_cache.json",
        "use_cache": True,
        "do_cws": False,
        "verbose": False,
        "nothrottle": False,
        "chunk_size": 0,
        "ksub": 0,
        "kwin": 0,
        "win_mode": 0,
        "min_hits": None,
        "cov_cws": 100,
        "refresh": 0,
        "simple": None,
        "full": False,
        "basis": None,
        "rounds": None,
        "out_dir": None,
        "out_dem": None,
        "out_stim": None,
        "morehelp": False,
        "version": False,
        "unrecognized": [],
    }

    def _has_later_input_file(start_idx: int) -> bool:
        for tok in argv[start_idx:]:
            if (
                not tok.startswith("-")
                and "=" not in tok
                and tok.endswith((".stim", ".dem", ".mtx", ".mmx"))
                and os.path.exists(tok)
            ):
                return True
        return False

    i = 0
    while i < len(argv):
        arg = argv[i]
        if not arg:
            i += 1
            continue

        if arg in ("--version", "-version", "version"):
            args["version"] = True
            i += 1
            continue

        if arg in ("--morehelp", "-morehelp", "--more-help", "-more-help", "morehelp"):
            args["morehelp"] = True
            i += 1
            continue

        if arg in ("-h", "--help", "help", "-?"):
            args["help"] = True
            i += 1
            continue

        if arg in ("--no-throttle", "-no-throttle", "nothrottle", "--nothrottle"):
            args["nothrottle"] = True
            i += 1
            continue

        if arg in ("-v", "--verbose", "verbose", "-verbose"):
            args["verbose"] = True
            i += 1
            continue

        if arg in ("--no-cache", "-no-cache", "nocache", "--nocache"):
            args["use_cache"] = False
            args["cache_file"] = None
            i += 1
            continue

        if arg in ("--trust-start", "-trust-start", "--trust_start", "-trust_start", "trust_start", "trust-start"):
            args["trust_start"] = True
            i += 1
            continue

        if arg in ("--cws", "-cws", "cws", "do_cws=1", "--do_cws"):
            args["do_cws"] = True
            i += 1
            continue

        if arg in (
            "--simple", "-simple", "--simplified", "-simplified",
            "simple", "simplified"
        ):
            args["simple"] = True
            args["full"] = False
            i += 1
            continue

        if arg in ("--full", "-full", "full"):
            args["full"] = True
            args["simple"] = False
            i += 1
            continue

        if arg in ("--rounds", "-rounds", "rounds"):
            if (
                i + 1 < len(argv)
                and "=" not in argv[i + 1]
                and argv[i + 1].lstrip("-").isdigit()
            ):
                args["rounds"] = int(argv[i + 1])
                i += 2
            else:
                args["rounds"] = 2
                i += 1
            continue

        if arg in ("--out-dem", "-out-dem", "--out_dem", "-out_dem", "out_dem"):
            next_tok = argv[i + 1] if i + 1 < len(argv) else None
            if (
                next_tok is not None
                and not next_tok.startswith("-")
                and "=" not in next_tok
                and not next_tok.endswith((".stim", ".mtx", ".mmx"))
                and not (
                    next_tok.endswith(".dem")
                    and os.path.exists(next_tok)
                    and args["fdem"] is None
                    and not _has_later_input_file(i + 2)
                )
            ):
                args["out_dem"] = next_tok
                i += 2
            else:
                args["out_dem"] = True
                i += 1
            continue

        if arg in ("--out-stim", "-out-stim", "--out_stim", "-out_stim", "out_stim"):
            next_tok = argv[i + 1] if i + 1 < len(argv) else None
            if (
                next_tok is not None
                and not next_tok.startswith("-")
                and "=" not in next_tok
                and not next_tok.endswith((".dem", ".mtx", ".mmx"))
                and not (
                    next_tok.endswith(".stim")
                    and os.path.exists(next_tok)
                    and args["fdem"] is None
                    and not _has_later_input_file(i + 2)
                )
            ):
                args["out_stim"] = next_tok
                i += 2
            else:
                args["out_stim"] = True
                i += 1
            continue

        key = None
        val = None
        if "=" in arg:
            key, val = arg.split("=", 1)
            if key.startswith("--"):
                key = key[2:]
            elif key.startswith("-"):
                key = key[1:]
        elif arg.startswith("--") or arg.startswith("-"):
            key = arg.lstrip("-")
            if i + 1 < len(argv) and not argv[i + 1].startswith("-") and "=" not in argv[i + 1]:
                val = argv[i + 1]
                i += 1
            else:
                val = "1"
        else:
            if os.path.exists(arg):
                if arg.endswith(".dem") or arg.endswith(".stim"):
                    args["fdem"] = arg
                elif arg.endswith(".mmx") or arg.endswith(".mtx"):
                    if args["finH"] is None:
                        args["finH"] = arg
                    elif args["finG"] is None and args["finL"] is None:
                        args["finG"] = arg
            else:
                args["unrecognized"].append(arg)
            i += 1
            continue

        if key:
            key_lower = key.lower()
            if key_lower in ("cache", "cache_file", "cachefile"):
                if val.lower() in ("0", "none", "false", "off", "no"):
                    args["use_cache"] = False
                    args["cache_file"] = None
                else:
                    args["use_cache"] = True
                    args["cache_file"] = val
            elif key_lower in ("fdem", "dem"):
                args["fdem"] = val
            elif key_lower == "finh":
                args["finH"] = val
            elif key_lower == "fing":
                args["finG"] = val
            elif key_lower == "finl":
                args["finL"] = val
            elif key_lower == "fin":
                args["fin"] = val
            elif key_lower in ("hx", "finhx"):
                args["Hx"] = val
            elif key_lower in ("hz", "finhz"):
                args["Hz"] = val
            elif key_lower in ("lx", "finlx"):
                args["Lx"] = val
            elif key_lower in ("lz", "finlz"):
                args["Lz"] = val
            elif key_lower == "finc":
                args["finC"] = val
            elif key_lower == "outc":
                args["outC"] = val
                args["do_cws"] = True
            elif key_lower in ("method", "m"):
                args["method"] = int(val)
            elif key_lower in ("dmin", "d_min"):
                args["dmin"] = int(val)
            elif key_lower in ("dmax", "d_max"):
                args["dmax"] = int(val)
            elif key_lower in ("wmin", "w_min"):
                args["wmin"] = int(val)
            elif key_lower in ("wmax", "w_max"):
                args["wmax"] = int(val)
            elif key_lower == "smax":
                args["smax"] = int(val)
            elif key_lower == "start":
                args["start"] = _normalize_start_list(val)
            elif key_lower in ("trust_start", "trust-start", "truststart"):
                args["trust_start"] = (
                    bool(int(val)) if val.isdigit() else (val.lower() not in ("0", "false", "no", "off"))
                )
            elif key_lower == "cbeg":
                args["cbeg"] = int(val)
            elif key_lower == "cend":
                args["cend"] = int(val)
            elif key_lower == "css":
                args["css"] = int(val)
            elif key_lower in ("dexp", "dest", "d_exp"):
                args["dexp"] = int(val)
            elif key_lower in ("steps", "num_steps", "nsteps"):
                args["steps"] = int(val)
            elif key_lower in ("threads", "num_threads", "t"):
                args["threads"] = int(val)
            elif key_lower in ("timeout", "time"):
                args["timeout"] = float(val)
            elif key_lower == "dw":
                args["dW"] = int(val)
            elif key_lower in ("maxc", "max_c"):
                args["maxC"] = int(val)
            elif key_lower == "pmin":
                args["pmin"] = float(val)
            elif key_lower == "noscan":
                args["noscan"] = int(val)
            elif key_lower == "classical":
                args["classical"] = int(val)
            elif key_lower == "seed":
                args["seed"] = int(val)
            elif key_lower in ("debug", "dbg"):
                args["debug"] |= parse_debug_value(val)  # OR-combined (the default is 0)
            elif key_lower == "solver":
                args["solver"] = val
            elif key_lower in ("verbose", "v"):
                args["verbose"] = bool(int(val)) if val.isdigit() else (val.lower() not in ("0", "false", "no", "off"))
            elif key_lower == "cws":
                args["do_cws"] = bool(int(val)) if val.isdigit() else (val.lower() not in ("0", "false", "no"))
            elif key_lower in ("nothrottle", "no_throttle"):
                args["nothrottle"] = (
                    bool(int(val)) if val.isdigit() else (val.lower() not in ("0", "false", "no", "off"))
                )
            elif key_lower in ("chunk_size", "chunksize", "batch", "chunk"):
                args["chunk_size"] = int(val)
            elif key_lower == "ksub":
                args["ksub"] = int(val)
            elif key_lower in ("kwin", "win"):
                args["kwin"] = int(val)
            elif key_lower in ("win_mode", "winmode"):
                args["win_mode"] = int(val)
            elif key_lower in ("min_hits", "minhits", "max_hits", "maxhits"):
                args["min_hits"] = int(val)
            elif key_lower in ("cov_cws", "covcws"):
                args["cov_cws"] = int(val)
            elif key_lower == "refresh":
                args["refresh"] = int(val)
            elif key_lower in ("simple", "simplified"):
                bval = (
                    bool(int(val)) if val.isdigit()
                    else (val.lower() not in ("0", "false", "no", "off"))
                )
                args["simple"] = bval
                if bval:
                    args["full"] = False
            elif key_lower == "full":
                bval = (
                    bool(int(val)) if val.isdigit()
                    else (val.lower() not in ("0", "false", "no", "off"))
                )
                args["full"] = bval
                if bval:
                    args["simple"] = False
            elif key_lower == "basis":
                if val.upper() not in ("X", "Z"):
                    raise ValueError(
                        f"Invalid basis '{val}'; must be 'X' or 'Z'."
                    )
                args["basis"] = val.upper()
            elif key_lower == "rounds":
                args["rounds"] = int(val)
            elif key_lower in ("out-dir", "out_dir", "outdir"):
                args["out_dir"] = val
            elif key_lower in ("out-dem", "out_dem", "outdem"):
                if val.lower() in ("0", "false", "no", "off", "none"):
                    args["out_dem"] = None
                elif val.lower() in ("1", "true", "yes", "auto", ""):
                    args["out_dem"] = True
                else:
                    args["out_dem"] = val
            elif key_lower in ("out-stim", "out_stim", "outstim"):
                if val.lower() in ("0", "false", "no", "off", "none"):
                    args["out_stim"] = None
                elif val.lower() in ("1", "true", "yes", "auto", ""):
                    args["out_stim"] = True
                else:
                    args["out_stim"] = val
            else:
                args["unrecognized"].append(arg)

        i += 1

    # Auto-infer classical mode when not explicitly set (matching src/util_io.c)
    if args["classical"] == -1:
        if (
            args["finG"] is not None or args["finL"] is not None
            or args["Hz"] is not None or args["Lz"] is not None
            or args["Lx"] is not None
            or args["fdem"] is not None or args["fin"] is not None
        ):
            args["classical"] = 0
        elif args["finH"] is not None or args["Hx"] is not None:
            args["classical"] = 1

    return args


def print_cli_short_help(file: Optional[Any] = None) -> None:
    text = f"""dist_m4ri.py (version {__version__}): Multithreaded distance calculator Python CLI
Usage: dist_m4ri.py [key=val | --flag val ...]

Allowed parameters:
  fdem, finH, finG, finL, fin, Hx, Hz, Lx, Lz, pmin, classical,
  --simple, --full, basis, --rounds, --out-dir, --out-dem, --out-stim,
  method, dmin, dmax, dexp (dest), steps, wmin, wmax, timeout, threads,
  nothrottle, chunk_size (batch), ksub, kwin (win), win_mode, min_hits,
  cov_cws, refresh, smax, noscan, start, cbeg, cend, finC, outC, maxC,
  dW, seed, debug, solver, cache, --no-cache, --verbose, --cws, --trust-start

Help options:
  -h, --help    : display help for commonly used parameters (fits 80 rows)
  --morehelp    : display full help for all available parameters
"""
    print(text, file=file)


def print_cli_help(file: Optional[Any] = None) -> None:
    help_text = f"""dist_m4ri.py (version {__version__}): Multithreaded distance calculator Python CLI
Usage: dist_m4ri.py [key=val | --flag val ...]

Input matrices & models:
  fdem=FILE             Detector Error Model (.dem) or Stim circuit (.stim) file
  finH=FILE             Parity check matrix H (classical) or Hx (CSS quantum) (.mmx/.mtx)
  finG=FILE, finL=FILE  Hz check matrix or Lx logical operator matrix (quantum CSS)
  fin=PREFIX            Base prefix for CSS matrices (loads ${{fin}}X.mtx, ${{fin}}Z.mtx, e.g. try -> tryX.mtx)
  Hx=FILE, Hz=FILE      CSS check matrices (alternative to finH/finG)
  Lx=FILE, Lz=FILE      CSS logical operators (optional, constructed if omitted)
  pmin=PROB             Minimum error probability threshold for DEM errors (default: 0.0)
  classical=0|1         1: classical code (Hx only), 0: quantum CSS (auto-detected)
  --simple / --full     For .stim files: keep primary-basis detectors only (--simple, default
                        for CSS/XZZX circuits) or keep all detectors (--full)
  basis=X|Z             Override primary memory basis ('X' or 'Z') for .stim circuits
  --rounds [N]          Set REPEAT block count in .stim circuit (default: 2; warns if no REPEAT)
  --out-dir DIR         Output directory for saved .dem/.stim files (default: input file dir)
  --out-dem [FILE]      Save constructed DEM (default: <basename>_simp.dem or _full.dem)
  --out-stim [FILE]     Save processed/noisy Stim circuit (default: <basename>_simp.stim or _full.stim)

Method and distance bounds:
  method=1|2|3          1=RW (upper bound), 2=CC (lower bound/exact), 3=Bracketing (default: 3)
  dmin=N                Certified lower bound, inclusive (CC starts from dmin) (default: 0)
  dmax=N                Known upper bound, inclusive (CC up to w=dmax-1 only) (default: 0)
  dexp=N                Expected distance estimate (alias: dest) (default: 0)

Search limits and stopping criteria:
  steps=N               Maximum RW steps / information sets (default: 100000)
  wmax=N                Maximum cluster weight to search in CC (0=until bound/timeout)
  wmin=N                Stop immediately if cw with weight <= wmin is found (default: 1)
  min_hits=N            Stop RW when lowest-weight cws are hit N times on average (default: 5)
  timeout=SEC           Execution timeout in seconds, 0 for infinite (default: 60.0)

Multithreading & execution:
  threads=N             Max worker threads to use (default: CPU cores; subject to
                        automatic throttling unless nothrottle=1)
  solver=NAME           Distance solver backend: 'dist_m4ri' (default) or 'codedistance'

Codeword collection & caching:
  --cws                 Collect and output non-trivial minimum-weight codewords
  outC=FILE             Save output codewords (for CSS, auto-suffixed _X.nz / _Z.nz)
  finC=FILE             Input initial candidate codewords (.nz file)
  cache=FILE            Persistent JSON cache file (default: tmp_dist_cache.json)
  --no-cache / nocache  Disable persistent JSON caching
  --verbose / -v        Explain bounds, steps, cache status, and miss probability estimates

Extra parameters (see --morehelp for details):
  smax=N (0)            Max syndrome weight for CC confinement profile (0 to disable)
  start=LIST (-1)       Expert: CC from listed columns only, e.g. start=0,48 (see --morehelp)
  --trust-start         Expert: accept start-list results as valid (copied to main cache record)
  noscan, cbeg, cend    Disabled in Python (ignored with a warning); use the dist_m4ri binary
  nothrottle=1 (0)      Disable automatic thread throttling (also --no-throttle)
  chunk_size=N (0)      RW batch chunk size (default: 0 for adaptive 25-500, alias: batch)
  ksub=N (0)            EXPERIMENTAL, do not use: RW subspace sketch dimension (0: full matrix)
  kwin=N (0)            RW localized window size (0: auto/hybrid for n>=500, alias: win)
  win_mode=0|1 (0)      RW window metric: 0=Tanner graph BFS, 1=index proximity
  cov_cws=N (100)       Number of lowest-weight cws used for the min_hits statistic
  refresh=N (0)         Periodic N basis refresh interval in RW steps (experimental ksub only)
  maxC=N (0)            Maximum number of codewords to collect (0: unlimited)
  dW=N (0)              Extra weight window above dmin to collect codewords
  seed=N (0)            Random number generator seed
  debug=N (0)           Debug bitmap, OR-combined (1: summary, 2: progress, 4: status, ...)

Help options:
  -h, --help            Display this help message (commonly used parameters)
  --morehelp            Display full help with all parameter descriptions
"""
    print(help_text, file=file)


def print_cli_morehelp(file: Optional[Any] = None) -> None:
    help_text = f"""dist_m4ri.py (version {__version__}): Multithreaded distance calculator Python CLI
Usage: dist_m4ri.py [key=val | --flag val ...]

Required input (at least one matrix/model specification):
  fdem=FILE             Detector Error Model (.dem) or Stim circuit (.stim) file.
                        Automatically constructs parity check H and logical L matrices.
  finH=FILE             Parity check matrix file in Matrix Market (.mmx / .mtx) format.
                        For classical codes, this is check matrix H. For CSS codes, Hx.
  finG=FILE             Generator / Hz matrix for quantum CSS codes in Matrix Market format.
  finL=FILE             Logical operator matrix Lx for quantum CSS codes in Matrix Market format.
                        Note: For a quantum CSS code, either finL (Lx) or finG (Hz) is required.
  fin=PREFIX            Base prefix for CSS matrices (loads ${{fin}}X.mtx and ${{fin}}Z.mtx, e.g. try -> tryX.mtx).
  Hx=FILE, Hz=FILE      Alternative syntax for specifying CSS check matrices Hx and Hz.
  Lx=FILE, Lz=FILE      Alternative syntax for specifying CSS logical operator matrices.
  pmin=PROB             Minimum error probability threshold for DEM parsing (default: 0.0).
                        Error mechanisms with probability < pmin are filtered out.
  classical=0|1         Code type override:
                        1: Classical linear code (Hx only; ignores/discards logicals).
                        0: Quantum CSS code (requires logicals or Hz).
                        Default: auto-detected (1 if only finH/Hx is given; 0 otherwise).
  --simple / simple=1   For .stim circuits: activate only primary-basis detectors by stripping
                        minority-basis detectors via ancilla and data basis tracking. Default
                        for CSS circuits (including locally rotated CSS such as XZZX codes).
  --full / full=1       For .stim circuits: retain all detectors (both X and Z sectors).
  basis=X|Z             Explicitly set primary memory basis ('X' or 'Z') for .stim circuits
                        (default: auto-detected from filename or initial/final resets/measurements).
  --rounds [N]          Construct DEM from a .stim circuit with N repetitions in its REPEAT block
                        (default: 2 when --rounds is given without a number; also rounds=N).
                        Issues a warning and continues with the original circuit if no REPEAT block.
  --out-dir DIR         Output directory for saving constructed DEM and/or processed Stim circuit
                        files (also out_dir=DIR; default: directory of the input .stim/.dem file).
  --out-dem [FILE]      Save constructed (flattened) DEM to FILE. If given without a filename (or
                        out_dem=1), automatically saves as <basename>_simp.dem or <basename>_full.dem.
  --out-stim [FILE]     Save processed (noisy, rounds-adjusted, detector-filtered) Stim circuit to
                        FILE. If given without a filename (or out_stim=1), automatically saves as
                        <basename>_simp.stim or <basename>_full.stim.

Calculation method:
  method=1|2|3          Calculation method (default: 3):
                        1: Random Window (RW) algorithm (upper bound on distance).
                           Repeatedly samples random information sets to find low-weight
                           codewords. Fast for discovering small errors.
                        2: Connected Cluster (CC) algorithm (lower bound / exact distance).
                           Exhaustive cluster search finding certified lower bound dmin or
                           exact distance if run to completion.
                        3: Bracketing mode (concurrent RW and CC).
                           Dynamically allocates worker threads between CC (lower bound) and
                           RW (upper bound) to converge on the exact distance rapidly, based on
                           the measured CC and RW times, the bounds, and the remaining timeout.

Distance bounds and guidance:
  dmin=N                Known certified lower bound on distance (default: 0).
                        CC search begins at weight w = max(1, dmin).
  dmax=N                Known upper bound on distance (default: 0).
                        CC ends once dmin = dmax (with outC, after the rounds w = dmax..dmax+dW);
                        RW ignores candidate codewords of weight >= dmax unless collecting codewords
                        or needed for min_hits.
  dexp=N                Expected code distance estimate (alias: dest) (default: 0).
                        Hint for method=3: as long as RW has found no codeword, CC rounds at
                        w > dexp run only on threads which cannot run RW (CC pauses if RW can
                        use all threads); CC resumes once RW finds a codeword or ends.

Search limits and stopping criteria:
  steps=N               Maximum number of RW steps / information sets across all threads
                        (default: 100000).
  wmax=N                Maximum cluster weight to analyze in CC (default: 0 = until bound/timeout;
                        with method=2, a known dmax is used as wmax).
  wmin=N                Minimum distance threshold (default: 1).
                        If a codeword of weight w <= wmin is discovered, search halts immediately
                        (not when collecting codewords with outC or maxC).
  timeout=SEC           Execution timeout in seconds (default: 60.0; set 0 for infinite).
                        In method=3, it guides the CC vs RW thread balance, and a CC round predicted
                        not to finish in time is not started while RW runs (in method=2, and with
                        steps=0, it is started anyway: it may still find a codeword of weight dmin).

Multithreading & throttling:
  threads=N             Maximum number of worker threads to allocate (default: min(CPU cores, 64)).
                        Subject to automatic thread throttling unless nothrottle=1 is specified:
                        - Small codes (<= 4 threads for n < 60, or for m*n < 2e4 with a small RW
                          workload m*n*steps < 5e7; otherwise <= 16 threads for n < 150 or
                          m*n < 1e5, if m*n*steps < 2e8; m: number of rows of H).
                        - Large memory matrices (dense RW working memory capped at ~1.5 GB).
                        - Small RW step counts (RW threads clamped to <= (steps + 9) / 10).
                        The last two limits apply to RW threads only: CC (method=2) is not
                        limited, and in method=3 CC rounds can use all threads.
  nothrottle=1          Disable automatic thread throttling (also --no-throttle).
                        Forces allocation of the exact number of threads requested.
  chunk_size=N          RW batch chunk size per thread (default: 0 = adaptive 25-500).
                        Alias: batch=N.
  ksub=N                EXPERIMENTAL, should not be used (a warning is printed): subspace sketch
                        dimension for RW (default: 0 = full matrix).
                        Precomputes N = ker(H) once; each thread echelonizes a compact
                        ksub x n sampled subspace in L1/L2 cache (e.g. ksub=32 or 64).
                        Automatically falls back to ksub=0 with a warning if m < nu = dim(ker H).
                        Each step reduces the span of ksub distinct rows of N, so that a few
                        codewords are found much more often than others, and min_hits may end
                        RW early.
  kwin=N                Localized column permutation window size W around a random seed
                        column j0 (default: 0 = automatic hybrid window/uniform for n>=500).
                        Alias: win=N.
  win_mode=0|1          Locality metric for kwin > 0: 0 = Tanner graph BFS neighbors
                        (default), 1 = contiguous column index window.
  min_hits=N            QDistRnd-style RW stopping criterion (default: 5, 0 = disabled).
                        Stops RW when the average number of times <n> that RW has found each
                        codeword reaches min_hits, both for the codewords of the minimum weight
                        found and for a representative set of cov_cws lowest-weight codewords
                        (heavier codewords are kept until enough lighter ones are found).  A
                        single codeword suffices.  A lighter codeword is missed with probability
                        about exp(-<n>), assuming that lighter codewords are found at least as
                        often.  With debug&1 (or --verbose), a warning is printed at the end if
                        the hit counts of codewords of the same weight are strongly non-uniform
                        (then exp(-<n>) is optimistic), together with an information-set estimate
                        (uniform random information sets) of the probability to miss a lighter
                        codeword and of the RW steps needed for 1%.  With --verbose, both
                        estimates are always shown for a result which is not certified.
  cov_cws=N             Size of the representative set of lowest-weight codewords for min_hits
                        (default: 100; 0 = only the minimum weight).  Without outC and maxC, at
                        most cov_cws minimum-weight codewords are tracked (all if cov_cws=0).
                        With fewer minimum-weight codewords, a smaller cov_cws may stop RW
                        earlier.
  refresh=N             Periodic adaptive basis refresh interval in RW steps, only with the
                        experimental ksub > 0 (default: 5000 when ksub > 0, 0 to disable).
                        Re-echelonizes N and substitutes heavier basis rows with discovered
                        min-weight codewords.

Connected Cluster (CC) search options:
  smax=N                Maximum syndrome weight for confinement profile (default: 0).
                        When smax > 0, tracks minimum syndrome weights for each error weight.
                        When smax=0, confinement is not computed.
  start=LIST            Expert option: comma-separated list of CC start columns (0-based), e.g.,
                        start=0,48,96 (default: -1 = all columns; start=N is a one-element list).
                        CC clusters are grown only from the listed columns, and they are NOT
                        limited to larger column indices.  Use, e.g., one column per block of a
                        quasi-cyclic code.  WARNING: dmin (and an exact result) is valid only if a
                        code symmetry maps every minimum-weight codeword to one containing a listed
                        column; dmax is always valid.  For DEM input, columns index the error
                        mechanisms in their order in the (flattened) DEM, after the pmin cutoff.
                        Results are cached under a separate key '<key>:start=a,b,c'; only dmax and
                        codewords are merged into the main cache record.  (Before version 0.10.0,
                        start=N limited CC to clusters started at column N using only larger
                        columns; this is now cbeg=N cend=N of the dist_m4ri binary.)
  --trust-start         Expert option (also trust_start=1): accept the start-list results as valid,
                        i.e., also merge dmin and rw_steps into the main cache record, and copy an
                        existing exact start-list record to the main record even if no calculation
                        is needed.  Only use it if the code symmetry assumption above holds.
  noscan, cbeg, cend    Disabled in the Python interface: accepted for backward compatibility, but
                        ignored with a warning to stderr.  These expert options of the dist_m4ri
                        binary do not give valid bounds on their own: noscan=1 runs CC only at
                        weight w=wmax (method=2); cbeg/cend limit CC start columns to a range for
                        split runs, where every column is eventually covered by one of the runs.
                        Run the dist_m4ri binary directly to use them (see its --morehelp).

Codeword collection and export:
  --cws                 Collect and display non-trivial minimum-weight codewords.
  outC=FILE             Export found codewords of weight up to min_w + dW (min_w: minimum weight
                        found) to file in .nz format.  In method=2/3, once the distance d is
                        known, CC rounds w = d..d+dW enumerate all such codewords (unless the
                        timeout is hit).
                        In CSS mode, automatically saves X-codewords to FILE_X.nz and
                        Z-codewords to FILE_Z.nz.
  finC=FILE             Import initial candidate codewords from file in .nz format.
                        In CSS mode, automatically resolves FILE_X.nz and FILE_Z.nz.
  maxC=N                Maximum number of codewords to collect (default: 0 = unlimited).
  dW=N                  Extra weight window above minimum distance to collect codewords
                        (w <= min_w + dW) (default: 0).

Distance caching (Python CLI):
  cache=FILE            Path to persistent JSON cache file (default: tmp_dist_cache.json).
  --no-cache / nocache  Disable reading and writing to the persistent JSON cache.

General options:
  solver=NAME           Distance calculation engine: 'dist_m4ri' (default) or 'codedistance'.
  --verbose / -v        Enable verbose output with detailed explanations of bounds,
                        timings, steps, and cache status (implies the debug bits 1, 2,
                        65536, and 131072).  For a RW upper bound which is not certified,
                        also estimates the probability that a lighter codeword was missed
                        (from the average hits <n>, and from uniform random information sets
                        with n, rank(H), and the RW steps), with the RW steps and time needed
                        for a 1% miss probability.
  seed=N                Random number generator seed (default: 0 = current time).
  debug=N               Debug bitmap (default: 0), a decimal or hexadecimal (0x...) integer;
                        multiple debug arguments are OR-combined.  Bits 1 to 32768 are passed
                        to the dist_m4ri binary, which prints diagnostic output to stderr (see
                        its --morehelp; lower bits are more informative):
                          1: summary (input matrices, stop reason, codeword statistics),
                          2: progress (run plan, finished CC rounds, new upper bounds),
                          4: periodic status (after 1, 2, 4, ... s, then every 60 s),
                          8: thread allocation and timing, 16: code parameters (ranks, k),
                          32: codeword supports, 64: arguments, 128: matrices,
                          256: codeword dump, 2048: no size cutoff.
                        Bits 65536 and above are used by the Python wrapper only:
                          65536: dist_m4ri command lines and run times,
                          131072: cache messages,
                          262144: Stim circuit to DEM conversion details,
                          524288: keep (and list) the temporary files.

Help options:
  -h, --help            Display summary help message (fits 80 rows).
  --morehelp            Display this full help message with all parameters.
"""
    print(help_text, file=file)


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    try:
        args = parse_cli_args(argv)
    except ValueError as e:
        sys.stderr.write(f"Error: {e}\n")
        return 1

    if args.get("version"):
        print(f"dist_m4ri.py version {__version__}")
        return 0

    if args.get("morehelp") or args.get("help"):
        compat_warn = check_binary_compatibility()
        if compat_warn:
            print(f"dist_m4ri.py: {compat_warn}", file=sys.stderr)
        if args.get("morehelp"):
            print_cli_morehelp()
        else:
            print_cli_help()
        return 0

    if args.get("unrecognized"):
        for u in args["unrecognized"]:
            print(f"dist_m4ri.py: unrecognized parameter '{u}'", file=sys.stderr)
        print_cli_short_help(file=sys.stderr)
        return 255

    if (
        not args.get("fdem") and not args.get("finH")
        and not args.get("Hx") and not args.get("Hz") and not args.get("fin")
    ):
        print("dist_m4ri.py: no input matrix or model specified", file=sys.stderr)
        print_cli_short_help(file=sys.stderr)
        return 255

    # When finC and outC are identical, empty or non-existent file is silently ignored (with a warning if verbose)
    args["finC"] = check_finc_outc(args["finC"], args["outC"], verbose=args["verbose"])

    cache_file = args["cache_file"] if args["use_cache"] else None
    if not args["use_cache"]:
        disable_distance_cache()

    try:
        if args["fdem"]:
            res = compute_dem_distance(
                dem=args["fdem"],
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                smax=args["smax"],
                start=args["start"],
                cbeg=args["cbeg"],
                cend=args["cend"],
                noscan=args["noscan"],
                dW=args["dW"],
                maxC=args["maxC"],
                pmin=args["pmin"],
                finC=args["finC"],
                outC=args["outC"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"],
                verbose=args["verbose"],
                nothrottle=args["nothrottle"],
                chunk_size=args["chunk_size"],
                ksub=args["ksub"],
                kwin=args["kwin"],
                win_mode=args["win_mode"],
                min_hits=args["min_hits"],
                cov_cws=args["cov_cws"],
                refresh=args["refresh"],
                simple=args["simple"],
                full=args["full"],
                basis=args["basis"],
                rounds=args["rounds"],
                out_dir=args["out_dir"],
                out_dem=args["out_dem"],
                out_stim=args["out_stim"],
                trust_start=args["trust_start"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, d_info, cws = res
                if args["outC"]:
                    _write_nzlist_file(args["outC"], cws)
            else:
                dist, d_info = res

            if args["verbose"]:
                print("=== DEM Distance Results ===")
                print(explain_bounds(
                    d_info, method=args["method"], label="DEM", stats=_last_run_stats
                ))
                print(f"  Summary bounds: {format_bounds_str(d_info)}")
            print(format_bounds_str(d_info))
            return 0

        if args.get("rounds") is not None:
            raise ValueError(
                "--rounds option is only supported for Stim circuit (.stim) inputs."
            )

        if args["Hx"] is not None or args["Hz"] is not None:
            res = compute_css_distance(
                Hx=args["Hx"],
                Hz=args["Hz"],
                Lx=args["Lx"],
                Lz=args["Lz"],
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                smax=args["smax"],
                start=args["start"],
                cbeg=args["cbeg"],
                cend=args["cend"],
                noscan=args["noscan"],
                dW=args["dW"],
                maxC=args["maxC"],
                finC=args["finC"],
                outC=args["outC"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"],
                verbose=args["verbose"],
                nothrottle=args["nothrottle"],
                chunk_size=args["chunk_size"],
                ksub=args["ksub"],
                kwin=args["kwin"],
                win_mode=args["win_mode"],
                min_hits=args["min_hits"],
                cov_cws=args["cov_cws"],
                refresh=args["refresh"],
                trust_start=args["trust_start"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, dx_info, dz_info, cws_x, cws_z = res
                if args["outC"]:
                    outC_X = _split_css_filename(args["outC"], "X")
                    outC_Z = _split_css_filename(args["outC"], "Z")
                    if cws_x:
                        _write_nzlist_file(outC_X, cws_x)
                    if cws_z:
                        _write_nzlist_file(outC_Z, cws_z)
            else:
                dist, dx_info, dz_info = res

            exact_tag = (
                " (exact)" if (
                    dx_info and dx_info[0] > 0 and dx_info[0] == dx_info[1]
                    and dz_info and dz_info[0] > 0 and dz_info[0] == dz_info[1]
                ) else ""
            )
            
            if args["verbose"]:
                print("=== CSS Quantum Code Distance Results ===")
                if dx_info:
                    print("--- X-Component Distance (dX) ---")
                    print(explain_bounds(
                        dx_info, method=args["method"], label="dX",
                        stats=_last_css_stats.get("X")
                    ))
                    print(f"  dX bounds: {format_bounds_str(dx_info)}")
                if dz_info:
                    print("--- Z-Component Distance (dZ) ---")
                    print(explain_bounds(
                        dz_info, method=args["method"], label="dZ",
                        stats=_last_css_stats.get("Z")
                    ))
                    print(f"  dZ bounds: {format_bounds_str(dz_info)}")
                print("--- Overall CSS Code Distance ---")
                print(f"  d = min(dX, dZ) = {dist}{exact_tag}")

            dx_str = format_bounds_str(dx_info) if dx_info else "none"
            dz_str = format_bounds_str(dz_info) if dz_info else "none"
            print(f"dX: {dx_str}  dZ: {dz_str}  (d = {dist}){exact_tag}")
            return 0

        # Handle fin prefix (e.g. fin=examples/try -> tryX.mtx and tryZ.mtx)
        finH = args["finH"]
        finG = args["finG"]
        if args["fin"]:
            if finH is None: finH = f"{args['fin']}X.mtx"
            if finG is None and args["finL"] is None: finG = f"{args['fin']}Z.mtx"

        # Quantum single-sided distance (finH with finG or finL, or classical=0)
        if finH and (finG is not None or args["finL"] is not None or args["classical"] == 0):
            res = compute_quantum_distance(
                H=finH,
                G=finG,
                L=args["finL"],
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                smax=args["smax"],
                start=args["start"],
                cbeg=args["cbeg"],
                cend=args["cend"],
                noscan=args["noscan"],
                dW=args["dW"],
                maxC=args["maxC"],
                finC=args["finC"],
                outC=args["outC"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                return_info=True,
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"],
                verbose=args["verbose"],
                nothrottle=args["nothrottle"],
                chunk_size=args["chunk_size"],
                ksub=args["ksub"],
                kwin=args["kwin"],
                win_mode=args["win_mode"],
                min_hits=args["min_hits"],
                cov_cws=args["cov_cws"],
                refresh=args["refresh"],
                trust_start=args["trust_start"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, d_info, cws = res
                if args["outC"]:
                    _write_nzlist_file(args["outC"], cws)
            else:
                dist, d_info = res

            if args["verbose"]:
                print("=== Quantum Code Distance Results (Single-Sided) ===")
                print(explain_bounds(
                    d_info, method=args["method"], label="Quantum", stats=_last_run_stats
                ))
                print(f"  Summary bounds: {format_bounds_str(d_info)}")
            print(format_bounds_str(d_info))
            return 0

        # Classical distance (finH only or classical=1)
        if finH:
            res = compute_classical_distance(
                H=finH,
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                smax=args["smax"],
                start=args["start"],
                cbeg=args["cbeg"],
                cend=args["cend"],
                noscan=args["noscan"],
                dW=args["dW"],
                maxC=args["maxC"],
                finC=args["finC"],
                outC=args["outC"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                return_info=True,
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"],
                verbose=args["verbose"],
                nothrottle=args["nothrottle"],
                chunk_size=args["chunk_size"],
                ksub=args["ksub"],
                kwin=args["kwin"],
                win_mode=args["win_mode"],
                min_hits=args["min_hits"],
                cov_cws=args["cov_cws"],
                refresh=args["refresh"],
                trust_start=args["trust_start"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, d_info, cws = res
                if args["outC"]:
                    _write_nzlist_file(args["outC"], cws)
            else:
                dist, d_info = res

            if args["verbose"]:
                print("=== Classical Code Distance Results ===")
                print(explain_bounds(
                    d_info, method=args["method"], label="Classical", stats=_last_run_stats
                ))
                print(f"  Summary bounds: {format_bounds_str(d_info)}")
            print(format_bounds_str(d_info))
            return 0
    except ValueError as e:
        sys.stderr.write(f"Error: {e}\n")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
