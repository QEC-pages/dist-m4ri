#!/usr/bin/env python3
"""Exclusive, restartable certification benchmark harness for `dist-m4ri` and future engines.

Implements the protocol of `scratch/bench_review/REVIEW.md` Sections 3.3 and 3.5 and writes
`dm4ri-bench-run/1`, `dm4ri-bench-env/1`, and `dm4ri-bench-plan/1` JSON files (see `SCHEMA.md`).

Key guarantees:
  - Input provenance (C1): verifies every tracked and cached input sha256 at session start and
    warms the OS page cache (m7).
  - Exclusive runs (C2): acquires `flock /tmp/qec-cert-bench.lock` (shared with qec-cert) around
    each timed run; `--no-lock` is provided strictly for isolated pinned testing.
  - Load gate & contamination (M4, m6): waits until <= 2 busy cores before 128-thread runs, uses
    `os.wait4` for child CPU time and max RSS, and measures `/proc/stat` busy/steal CPU-seconds
    to flag runs where external CPU time exceeds 2% of `threads * wall_sec`.
  - Budgets, seeds, placement, kill (M1, M7, M8, m9, m10): per-tier budgets 60/300/1800/3600 s,
    seeds 5/5/3/1 (always >= 1), p=128 + p=64 on T1/T2 (p=1 pinned to CPU 1, not CPU 0), fixed-
    seed random permutation per tier, periodic canary run, SIGTERM -> SIGKILL at
    T + max(30 s, 0.05 T).
  - Independent verification (C4): invokes `verify.py` after every run outside the timed region.
  - Restartable & atomic: `run_id = sha1(canonical_json(key))`; finished runs are skipped on
    restart (`--retry` retries harness errors).

Usage:
  python3 tools/bench/harness.py --manifest benchmark/manifest/certify_v1.json --session <dir>
  python3 tools/bench/harness.py --session <dir> --plan-only
  python3 tools/bench/harness.py --session <dir> --tiers T0 --threads 2 --seeds 1 \\
      --taskset 126-127 --nice 19 --no-lock --no-load-gate --canary-every 0
"""
from __future__ import annotations

import argparse
import datetime
import fcntl
import hashlib
import json
import os
import platform
import random
import re
import resource
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

HARNESS_VERSION = "1.0.0"
HARNESS_MAJOR = 1
DEFAULT_LOCK_FILE = os.environ.get("QEC_CERT_BENCH_LOCK", "/tmp/qec-cert-bench.lock")

BENCH_DIR = Path(__file__).resolve().parent
REPO_ROOT = BENCH_DIR.parent.parent
MANIFEST_DIR = REPO_ROOT / "benchmark" / "manifest"

for _p in (str(BENCH_DIR), str(MANIFEST_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import make_inputs as mi  # noqa: E402
import verify as vf  # noqa: E402
from adapters import Adapter, build_adapter, git_info  # noqa: E402
from adapters.dist_m4ri_py import find_qec_cert_lib, find_qec_cert_root  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers: canonical JSON, atomic IO, /proc/stat, env.json
# ---------------------------------------------------------------------------

def canonical_json(obj: Any) -> str:
    """Deterministic compact JSON (sorted keys, no extra whitespace)."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def compute_run_id(key: Dict[str, Any]) -> str:
    """Deterministic sha1 hex id for a planned run key."""
    return hashlib.sha1(canonical_json(key).encode("utf-8")).hexdigest()


def utc_now_iso() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def atomic_write_json(path: Path, obj: Any, *, indent: Optional[int] = 2) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    text = (json.dumps(obj, indent=indent) if indent is not None else json.dumps(obj)) + "\n"
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def append_jsonl(path: Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(obj, separators=(",", ":")) + "\n"
    with open(path, "a", encoding="utf-8") as f:
        f.write(line)
        f.flush()
        os.fsync(f.fileno())


def read_proc_stat_jiffies() -> Tuple[int, int, int]:
    """Return (busy_jiffies, steal_jiffies, total_jiffies) from aggregate cpu line of /proc/stat.

    Linux /proc/stat fields: user, nice, system, idle, iowait, irq, softirq, steal[, guest, guest_nice].
    `guest` and `guest_nice` are already accounted inside `user` and `nice`, so we sum only fields 0..7.
    """
    with open("/proc/stat", "r", encoding="utf-8") as f:
        parts = f.readline().split()[1:]
    vals = [int(x) for x in parts[:8]]
    while len(vals) < 8:
        vals.append(0)
    user, nice, system, idle, iowait, irq, softirq, steal = vals
    busy = user + nice + system + irq + softirq + steal
    total = busy + idle + iowait
    return busy, steal, total


def measure_busy_cores(window_s: float = 2.0) -> float:
    """Measure machine-wide busy cores over `window_s` seconds."""
    b0, _s0, t0 = read_proc_stat_jiffies()
    time.sleep(max(0.1, window_s))
    b1, _s1, t1 = read_proc_stat_jiffies()
    ncpu = os.cpu_count() or 128
    dt = max(1, t1 - t0)
    return round(ncpu * max(0, b1 - b0) / dt, 3)


def default_taskset_for_threads(threads: int) -> Optional[str]:
    """Thread placement (REVIEW.md M7): p=1 avoids CPU 0; p=64 uses physical cores 0-63; p=128 uses 0-127."""
    ncpu = os.cpu_count() or 128
    if threads <= 0:
        return None
    if threads == 1:
        return "1" if ncpu > 1 else "0"
    hi = min(threads, ncpu) - 1
    return f"0-{hi}"


def _run_cmd_text(cmd: Sequence[str], timeout: float = 5.0) -> Optional[str]:
    try:
        p = subprocess.run(list(cmd), capture_output=True, text=True, timeout=timeout)
        if p.returncode == 0:
            return p.stdout.strip()
    except Exception:
        pass
    return None


def _parse_lscpu(raw: str) -> Dict[str, Any]:
    kv: Dict[str, str] = {}
    numa: Dict[str, str] = {}
    caches: Dict[str, str] = {}
    for line in raw.splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        k, v = k.strip(), v.strip()
        kv[k] = v
        if k.startswith("NUMA node") and k.endswith("CPU(s)"):
            numa[k] = v
        elif k.endswith("cache"):
            caches[k] = v
    flags = kv.get("Flags", "")
    flags_sha1 = hashlib.sha1(flags.encode("utf-8")).hexdigest() if flags else None
    return {
        "model": kv.get("Model name"),
        "architecture": kv.get("Architecture"),
        "cpus": int(kv["CPU(s)"]) if kv.get("CPU(s)", "").isdigit() else os.cpu_count(),
        "sockets": int(kv["Socket(s)"]) if kv.get("Socket(s)", "").isdigit() else None,
        "cores_per_socket": int(kv["Core(s) per socket"]) if kv.get("Core(s) per socket", "").isdigit() else None,
        "threads_per_core": int(kv["Thread(s) per core"]) if kv.get("Thread(s) per core", "").isdigit() else None,
        "numa": numa,
        "caches": caches,
        "hypervisor": kv.get("Hypervisor vendor"),
        "flags_sha1": flags_sha1,
        "lscpu": raw,
    }


def _parse_makefile_flags(repo_root: Path) -> Dict[str, Any]:
    mk = repo_root / "src" / "makefile"
    flags: Dict[str, str] = {}
    if mk.is_file():
        for line in mk.read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            for key in ("CC", "OPT", "WARN", "EXTRA", "CFLAGS"):
                if re.match(rf"^{key}\s*=", s):
                    flags[key] = s.split("=", 1)[1].strip()
    cc_bin = flags.get("CC", "gcc")
    cc_ver = _run_cmd_text([cc_bin, "--version"])
    return {
        "cc": cc_ver.splitlines()[0] if cc_ver else cc_bin,
        "flags": flags,
    }


def _pkg_version(name: str) -> Optional[str]:
    try:
        import importlib.metadata

        return importlib.metadata.version(name)
    except Exception:
        try:
            mod = __import__(name)
            return str(getattr(mod, "__version__", "unknown"))
        except Exception:
            return None


def _read_first_line(path: Path) -> Optional[str]:
    try:
        if path.is_file():
            return path.read_text(encoding="utf-8").splitlines()[0].strip()
    except Exception:
        pass
    return None


def collect_env_json(
    *,
    repo_root: Path,
    binary_path: Path,
    manifest: Dict[str, Any],
    manifest_sha: str,
    inputs_summary: Dict[str, Any],
    qec_cert_root: Optional[Path],
) -> Tuple[Dict[str, Any], str]:
    """Collect `dm4ri-bench-env/1` environment provenance and its sha1 id."""
    lscpu_raw = _run_cmd_text(["lscpu"]) or ""
    mem_total_kb: Optional[int] = None
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                mem_total_kb = int(line.split()[1])
                break
    except Exception:
        pass

    os_pretty = None
    try:
        for line in Path("/etc/os-release").read_text(encoding="utf-8").splitlines():
            if line.startswith("PRETTY_NAME="):
                os_pretty = line.split("=", 1)[1].strip().strip('"')
                break
    except Exception:
        pass

    autogroup = _read_first_line(Path("/proc/sys/kernel/sched_autogroup_enabled"))
    vm_product = _read_first_line(Path("/sys/class/dmi/id/product_name"))

    dm_commit, dm_dirty = git_info(repo_root)
    dm_bin_ver = None
    if binary_path.is_file():
        vtxt = _run_cmd_text([str(binary_path), "--version"])
        if vtxt:
            m = re.search(r"([0-9]+\.[0-9]+\.[0-9]+[^\s]*)", vtxt)
            dm_bin_ver = m.group(1) if m else vtxt
    py_path = repo_root / "dist_m4ri.py"
    py_ver = None
    if py_path.is_file():
        for line in py_path.read_text(encoding="utf-8").splitlines()[:80]:
            m = re.match(r'^__version__\s*=\s*["\']([^"\']+)["\']', line.strip())
            if m:
                py_ver = m.group(1)
                break

    q_root = find_qec_cert_root(repo_root, explicit=qec_cert_root)
    q_commit, q_dirty = git_info(q_root) if q_root else (None, None)
    q_lib = find_qec_cert_lib(q_root)
    q_ver = None
    if q_root is not None:
        for rel in ("qec_cert/__init__.py", "pyproject.toml"):
            p = q_root / rel
            if p.is_file():
                for line in p.read_text(encoding="utf-8").splitlines()[:60]:
                    m = re.match(r'^(?:__version__|version)\s*=\s*["\']([^"\']+)["\']', line.strip())
                    if m:
                        q_ver = m.group(1)
                        break
            if q_ver:
                break

    harness_files: Dict[str, str] = {}
    for p in sorted(BENCH_DIR.rglob("*.py")):
        rel = str(p.relative_to(REPO_ROOT))
        harness_files[rel] = mi.sha256_file(p)
    mk_inp = MANIFEST_DIR / "make_inputs.py"
    if mk_inp.is_file():
        harness_files[str(mk_inp.relative_to(REPO_ROOT))] = mi.sha256_file(mk_inp)

    env_obj: Dict[str, Any] = {
        "schema": "dm4ri-bench-env/1",
        "created": utc_now_iso(),
        "host": socket.gethostname(),
        "kernel": platform.release(),
        "os": os_pretty or platform.platform(),
        "vm": vm_product,
        "sched": {"autogroup": int(autogroup) if (autogroup and autogroup.isdigit()) else autogroup},
        "cpu": _parse_lscpu(lscpu_raw),
        "memory": {"total_kb": mem_total_kb},
        "compiler": _parse_makefile_flags(repo_root),
        "m4ri": _run_cmd_text(["pkg-config", "--modversion", "m4ri"]),
        "python": {
            "version": sys.version.split()[0],
            "executable": sys.executable,
            "stim": _pkg_version("stim"),
            "numpy": _pkg_version("numpy"),
            "scipy": _pkg_version("scipy"),
            "matplotlib": _pkg_version("matplotlib"),
        },
        "qec_cert": {
            "path": str(q_root) if q_root else None,
            "version": q_ver,
            "commit": q_commit,
            "dirty": q_dirty,
            "lib": str(q_lib) if q_lib else None,
            "lib_sha256": mi.sha256_file(q_lib) if q_lib else None,
        },
        "dist_m4ri": {
            "repo": str(repo_root),
            "version": dm_bin_ver,
            "commit": dm_commit,
            "dirty": dm_dirty,
            "binary": str(binary_path),
            "sha256": mi.sha256_file(binary_path) if binary_path.is_file() else None,
            "py": str(py_path),
            "py_version": py_ver,
            "py_sha256": mi.sha256_file(py_path) if py_path.is_file() else None,
        },
        "harness": {"version": HARNESS_VERSION, "files": harness_files},
        "manifest": {"id": manifest["id"], "version": manifest["version"], "sha256": manifest_sha},
        "inputs": inputs_summary,
    }
    env_id = hashlib.sha1(canonical_json(env_obj).encode("utf-8")).hexdigest()
    return env_obj, env_id


# ---------------------------------------------------------------------------
# Process execution with wait4, flock, load gate, kill at T + max(30, 0.05 T)
# ---------------------------------------------------------------------------

def kill_after_seconds(budget_s: float) -> float:
    """Kill deadline per REVIEW.md 3.3: T + max(30 s, 0.05 T)."""
    return float(budget_s) + max(30.0, 0.05 * float(budget_s))


def run_process_wait4(
    argv: Sequence[str],
    *,
    cwd: str,
    env_extra: Optional[Dict[str, str]],
    out_path: Path,
    err_path: Path,
    kill_after_s: float,
    taskset_cpus: Optional[str] = None,
    nice_val: Optional[int] = None,
    sigterm_grace_s: float = 10.0,
) -> Dict[str, Any]:
    """Run one subprocess in its own session, wait with `os.wait4`, and measure CPU/RSS/contamination."""
    cmd: List[str] = []
    if taskset_cpus:
        cmd.extend(["taskset", "-c", taskset_cpus])
    if nice_val is not None:
        cmd.extend(["nice", "-n", str(nice_val)])
    cmd.extend(argv)

    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)

    clk_tck = float(os.sysconf("SC_CLK_TCK"))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    err_path.parent.mkdir(parents=True, exist_ok=True)

    started_at = utc_now_iso()
    busy0, steal0, _tot0 = read_proc_stat_jiffies()
    ru_self0 = resource.getrusage(resource.RUSAGE_SELF)
    t0 = time.monotonic()

    killed = False
    term_sent_at: Optional[float] = None
    kill_sent = False
    status_raw = 0
    rusage = None

    with open(out_path, "wb") as out_f, open(err_path, "wb") as err_f:
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            env=env,
            stdout=out_f,
            stderr=err_f,
            start_new_session=True,
        )
        pid = proc.pid
        sleep_s = 0.01
        while True:
            wpid, st, ru = os.wait4(pid, os.WNOHANG)
            if wpid == pid:
                status_raw = st
                rusage = ru
                break
            elapsed = time.monotonic() - t0
            if not killed and elapsed >= kill_after_s:
                killed = True
                term_sent_at = time.monotonic()
                try:
                    os.killpg(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            elif killed and not kill_sent and term_sent_at is not None:
                if (time.monotonic() - term_sent_at) >= sigterm_grace_s:
                    kill_sent = True
                    try:
                        os.killpg(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
            time.sleep(sleep_s)
            sleep_s = min(0.2, sleep_s * 1.25)

    wall_sec = time.monotonic() - t0
    ended_at = utc_now_iso()
    busy1, steal1, _tot1 = read_proc_stat_jiffies()
    ru_self1 = resource.getrusage(resource.RUSAGE_SELF)

    # Ensure Python's Popen object knows the child was reaped
    proc.returncode = os.waitstatus_to_exitcode(status_raw)
    if os.WIFEXITED(status_raw):
        rc: Optional[int] = os.WEXITSTATUS(status_raw)
        sig: Optional[int] = None
    elif os.WIFSIGNALED(status_raw):
        rc = None
        sig = os.WTERMSIG(status_raw)
    else:
        rc = proc.returncode
        sig = None

    cpu_user_s = float(rusage.ru_utime) if rusage is not None else 0.0
    cpu_sys_s = float(rusage.ru_stime) if rusage is not None else 0.0
    max_rss_kb = int(rusage.ru_maxrss) if rusage is not None else 0

    busy_cpu_s = max(0.0, (busy1 - busy0) / clk_tck)
    steal_s = max(0.0, (steal1 - steal0) / clk_tck)
    child_cpu_s = cpu_user_s + cpu_sys_s
    harness_cpu_s = max(
        0.0, (ru_self1.ru_utime + ru_self1.ru_stime) - (ru_self0.ru_utime + ru_self0.ru_stime)
    )
    ext_cpu_s = max(0.0, busy_cpu_s - child_cpu_s - harness_cpu_s)

    return {
        "cmd": cmd,
        "rc": rc,
        "signal": sig,
        "killed": killed,
        "wall_sec": round(wall_sec, 6),
        "cpu_user_s": round(cpu_user_s, 6),
        "cpu_sys_s": round(cpu_sys_s, 6),
        "max_rss_kb": max_rss_kb,
        "busy_cpu_s": round(busy_cpu_s, 4),
        "child_cpu_s": round(child_cpu_s, 4),
        "harness_cpu_s": round(harness_cpu_s, 4),
        "ext_cpu_s": round(ext_cpu_s, 4),
        "steal_s": round(steal_s, 4),
        "started_at": started_at,
        "ended_at": ended_at,
    }


def wait_for_load_gate(
    *,
    max_busy: float,
    window_s: float,
    poll_s: float,
    max_wait_s: Optional[float],
    log_fn: Any,
) -> Dict[str, Any]:
    """Wait until external machine busy cores <= max_busy (two consecutive readings after pausing)."""
    t0 = time.monotonic()
    busy = measure_busy_cores(window_s)
    if busy <= max_busy:
        return {"busy_cores": busy, "max_busy": max_busy, "waited_s": 0.0, "passed": True}
    log_fn(f"LOAD GATE: busy_cores={busy:.2f} > {max_busy:.2f}; waiting")
    consec_ok = 0
    while True:
        waited = time.monotonic() - t0
        if max_wait_s is not None and waited >= max_wait_s:
            log_fn(f"LOAD GATE: exceeded max wait ({waited:.1f}s); busy_cores={busy:.2f}")
            return {"busy_cores": busy, "max_busy": max_busy, "waited_s": round(waited, 2), "passed": False}
        time.sleep(poll_s if consec_ok == 0 else min(10.0, poll_s))
        busy = measure_busy_cores(window_s)
        if busy <= max_busy:
            consec_ok += 1
            if consec_ok >= 2:
                waited = time.monotonic() - t0
                log_fn(f"LOAD GATE: cleared after {waited:.1f}s (busy_cores={busy:.2f})")
                return {"busy_cores": busy, "max_busy": max_busy, "waited_s": round(waited, 2), "passed": True}
        else:
            consec_ok = 0


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

def build_plan(
    manifest: Dict[str, Any],
    adapters: Sequence[Adapter],
    *,
    track: str = "certify",
    tiers: Optional[Sequence[str]] = None,
    instance_ids: Optional[Sequence[str]] = None,
    only_small: bool = False,
    threads_override: Optional[Sequence[int]] = None,
    seeds_override: Optional[Sequence[int]] = None,
    budget_override: Optional[float] = None,
    order_seed: int = 20261124,
) -> Dict[str, Any]:
    """Build deterministic randomized run plan (`dm4ri-bench-plan/1`)."""
    rng = random.Random(order_seed)
    tier_order = ["T0", "T1", "T2", "T3"]
    if tiers:
        tier_order = [t for t in tier_order if t in set(tiers)]

    if track == "scaling":
        league = "scaling"
        target_ids = list(instance_ids) if instance_ids else list(manifest.get("scaling_instances") or [])
        selected = mi.select_instances(manifest, ids=target_ids, only_small=only_small)
        p_list = list(threads_override) if threads_override else [128, 64, 32, 16, 8, 4, 2, 1]
        planned_runs: List[Dict[str, Any]] = []
        for inst in selected:
            tier_cfg = manifest["tiers"][inst["tier"]]
            budget = float(budget_override) if budget_override is not None else float(tier_cfg["budget_s"])
            for ad in adapters:
                if not ad.supports(inst):
                    continue
                seeds = list(seeds_override) if seeds_override else [1, 2, 3]
                for p in p_list:
                    for seed in seeds:
                        if seed < 1:
                            raise ValueError(f"Seed must be >= 1 (got {seed})")
                        key = {
                            "manifest": f"{manifest['id']}/{manifest['version']}",
                            "instance": inst["id"],
                            "engine": ad.name,
                            "config": ad.config,
                            "league": league,
                            "threads": int(p),
                            "budget_s": budget,
                            "seed": int(seed),
                            "canary": None,
                            "harness_major": HARNESS_MAJOR,
                        }
                        planned_runs.append({
                            "run_id": compute_run_id(key),
                            "key": key,
                            "tier": inst["tier"],
                            "instance": inst["id"],
                            "engine": ad.name,
                            "config": ad.config,
                            "league": league,
                            "threads": int(p),
                            "budget_s": budget,
                            "seed": int(seed),
                        })
        rng.shuffle(planned_runs)
    else:
        league = "default"
        selected = mi.select_instances(manifest, tiers=tier_order, ids=instance_ids, only_small=only_small)
        planned_runs = []
        for tier in tier_order:
            tier_insts = [inst for inst in selected if inst["tier"] == tier]
            if not tier_insts:
                continue
            tier_cfg = manifest["tiers"][tier]
            budget = float(budget_override) if budget_override is not None else float(tier_cfg["budget_s"])
            seeds = [int(s) for s in (seeds_override if seeds_override else tier_cfg["seeds"])]
            if any(s < 1 for s in seeds):
                raise ValueError(f"Harness rule m9: all seeds must be >= 1 (got {seeds})")

            tier_runs: List[Dict[str, Any]] = []
            for inst in tier_insts:
                for ad in adapters:
                    if not ad.supports(inst):
                        continue
                    if threads_override is not None:
                        p_list = [int(x) for x in threads_override]
                    elif ad.name == "dist_m4ri" and ad.config == "default":
                        # 128 on all tiers + 64 on T1/T2 (from manifest tiers table)
                        p_list = [int(x) for x in tier_cfg["threads"]]
                    else:
                        p_list = [int(tier_cfg["threads"][0])]
                    for p in p_list:
                        for seed in seeds:
                            key = {
                                "manifest": f"{manifest['id']}/{manifest['version']}",
                                "instance": inst["id"],
                                "engine": ad.name,
                                "config": ad.config,
                                "league": league,
                                "threads": int(p),
                                "budget_s": budget,
                                "seed": int(seed),
                                "canary": None,
                                "harness_major": HARNESS_MAJOR,
                            }
                            tier_runs.append({
                                "run_id": compute_run_id(key),
                                "key": key,
                                "tier": tier,
                                "instance": inst["id"],
                                "engine": ad.name,
                                "config": ad.config,
                                "league": league,
                                "threads": int(p),
                                "budget_s": budget,
                                "seed": int(seed),
                            })
            rng.shuffle(tier_runs)
            planned_runs.extend(tier_runs)

    return {
        "schema": "dm4ri-bench-plan/1",
        "created": utc_now_iso(),
        "order_seed": order_seed,
        "track": track,
        "manifest": {"id": manifest["id"], "version": manifest["version"]},
        "options": {
            "tiers": tier_order,
            "only_small": only_small,
            "engines": [f"{a.name}:{a.config}" for a in adapters],
            "threads_override": list(threads_override) if threads_override else None,
            "seeds_override": list(seeds_override) if seeds_override else None,
            "budget_override": budget_override,
        },
        "total_runs": len(planned_runs),
        "runs": planned_runs,
    }


# ---------------------------------------------------------------------------
# Session execution
# ---------------------------------------------------------------------------

def load_existing_records(results_path: Path) -> Dict[str, Dict[str, Any]]:
    """Load results.jsonl into {run_id: last_record}."""
    done: Dict[str, Dict[str, Any]] = {}
    if not results_path.is_file():
        return done
    for line in results_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
            rid = rec.get("run_id")
            if rid:
                done[str(rid)] = rec
        except ValueError:
            pass
    return done


def should_skip_run(existing: Optional[Dict[str, Any]], *, retry: bool, retry_all_errors: bool) -> bool:
    if existing is None:
        return False
    res = existing.get("result") or {}
    if res.get("status") != "error":
        return True
    err_kind = (res.get("error") or {}).get("kind")
    if retry_all_errors:
        return False
    if retry and err_kind == "harness":
        return False
    return True


def execute_single_run(
    *,
    item: Dict[str, Any],
    inst: Dict[str, Any],
    adapter: Adapter,
    verifier: vf.VerifierContext,
    session_dir: Path,
    session_id: str,
    manifest: Dict[str, Any],
    manifest_sha: str,
    env_id: str,
    lock_file: Optional[str],
    use_load_gate: bool,
    max_busy_cores: float,
    load_window_s: float,
    load_poll_s: float,
    load_max_wait_s: Optional[float],
    taskset_override: Optional[str],
    nice_val: Optional[int],
    log_fn: Any,
) -> Dict[str, Any]:
    """Execute one run (lock -> load gate -> prepare -> wait4 subprocess -> parse -> verify -> record)."""
    run_id = item["run_id"]
    key = item["key"]
    threads = int(item["threads"])
    budget_s = float(item["budget_s"])
    seed = int(item["seed"])
    canary = key.get("canary")
    kill_s = kill_after_seconds(budget_s)
    run_dir = session_dir / "runs" / run_id
    out_path = run_dir / "out.txt"
    err_path = run_dir / "err.txt"
    taskset_cpus = taskset_override if taskset_override is not None else default_taskset_for_threads(threads)

    lock_fd = None
    lock_wait_s = 0.0
    load_gate_info: Dict[str, Any] = {
        "busy_cores": None,
        "max_busy": max_busy_cores if use_load_gate else None,
        "waited_s": 0.0,
        "passed": True,
    }

    try:
        if lock_file:
            t_lk0 = time.monotonic()
            lock_fd = os.open(lock_file, os.O_CREAT | os.O_RDWR, 0o666)
            fcntl.flock(lock_fd, fcntl.LOCK_EX)
            lock_wait_s = round(time.monotonic() - t_lk0, 4)

        if use_load_gate:
            load_gate_info = wait_for_load_gate(
                max_busy=max_busy_cores,
                window_s=load_window_s,
                poll_s=load_poll_s,
                max_wait_s=load_max_wait_s,
                log_fn=log_fn,
            )

        prep = adapter.prepare(inst)
        prep_s = float(prep.get("prep_s", 0.0))
        cmd_spec = adapter.command(inst, threads, budget_s, seed, {"prep": prep})

        proc_info = run_process_wait4(
            cmd_spec["argv"],
            cwd=cmd_spec.get("cwd") or str(REPO_ROOT),
            env_extra=cmd_spec.get("env"),
            out_path=out_path,
            err_path=err_path,
            kill_after_s=kill_s,
            taskset_cpus=taskset_cpus,
            nice_val=nice_val,
        )
    except Exception as exc:
        # Harness-level error: record it so `--retry` can pick it up cleanly
        now = utc_now_iso()
        record = _assemble_error_record(
            item=item,
            inst=inst,
            adapter=adapter,
            session_id=session_id,
            manifest=manifest,
            manifest_sha=manifest_sha,
            env_id=env_id,
            taskset_cpus=taskset_cpus,
            nice_val=nice_val,
            lock_file=lock_file,
            lock_wait_s=lock_wait_s,
            load_gate_info=load_gate_info,
            kill_s=kill_s,
            run_dir=run_dir,
            session_dir=session_dir,
            err_msg=f"{type(exc).__name__}: {exc}",
            timestamp=now,
        )
        atomic_write_json(run_dir / "run.json", record)
        append_jsonl(session_dir / "results.jsonl", record)
        return record
    finally:
        if lock_fd is not None:
            try:
                fcntl.flock(lock_fd, fcntl.LOCK_UN)
            finally:
                os.close(lock_fd)

    # Parse and verify OUTSIDE the lock and timed region
    stdout_text = out_path.read_text(encoding="utf-8", errors="replace") if out_path.is_file() else ""
    stderr_text = err_path.read_text(encoding="utf-8", errors="replace") if err_path.is_file() else ""
    parsed = adapter.parse(
        proc_info["rc"],
        stdout_text,
        stderr_text,
        killed=proc_info["killed"],
        prep=prep,
    )

    wall_sec = proc_info["wall_sec"]
    contaminated = bool(proc_info["ext_cpu_s"] > 0.02 * threads * max(wall_sec, 0.001))
    exclusive = bool(lock_file is not None and (not use_load_gate or load_gate_info.get("passed")))

    result_block = dict(parsed)
    result_block.update(
        wall_sec=wall_sec,
        cap_sec=budget_s,
        prep_s=round(prep_s, 6),
        cpu_user_s=proc_info["cpu_user_s"],
        cpu_sys_s=proc_info["cpu_sys_s"],
        max_rss_kb=proc_info["max_rss_kb"],
        rc=proc_info["rc"],
        signal=proc_info["signal"],
        killed=proc_info["killed"],
    )

    record: Dict[str, Any] = {
        "schema": "dm4ri-bench-run/1",
        "run_id": run_id,
        "key": key,
        "session": session_id,
        "harness": {"version": HARNESS_VERSION},
        "manifest": {"id": manifest["id"], "version": manifest["version"], "sha256": manifest_sha},
        "instance": {
            "id": inst["id"],
            "kind": inst["kind"],
            "tier": inst["tier"],
            "d_ref": inst.get("d_ref"),
            "bracket_ref": inst.get("bracket_ref"),
            "compare3_id": inst.get("compare3_id"),
        },
        "engine": adapter.describe(),
        "options": {
            "league": key["league"],
            "threads": threads,
            "budget_s": budget_s,
            "seed": seed,
            "argv_extra": cmd_spec.get("argv_extra") or [],
            "hints": {},
            "canary": canary,
            "kill_after_s": round(kill_s, 3),
            "cpus": taskset_cpus,
        },
        "placement": {
            "taskset": taskset_cpus,
            "nice": nice_val,
            "exclusive": exclusive,
            "lock": lock_file,
            "lock_wait_sec": lock_wait_s,
            "load_gate": load_gate_info,
        },
        "result": result_block,
        "verify": {},
        "env": {
            "env_id": env_id,
            "busy_cpu_s": proc_info["busy_cpu_s"],
            "child_cpu_s": proc_info["child_cpu_s"],
            "harness_cpu_s": proc_info["harness_cpu_s"],
            "ext_cpu_s": proc_info["ext_cpu_s"],
            "steal_s": proc_info["steal_s"],
            "contaminated": contaminated,
            "started_at": proc_info["started_at"],
            "ended_at": proc_info["ended_at"],
        },
        "raw": {
            "dir": str(run_dir.relative_to(session_dir)),
            "stdout": str(out_path.relative_to(session_dir)),
            "stderr": str(err_path.relative_to(session_dir)),
            "argv": proc_info["cmd"],
        },
    }

    verifier.verify_record(record, inst, dem_path_hint=prep.get("dem_path"))
    atomic_write_json(run_dir / "run.json", record)
    append_jsonl(session_dir / "results.jsonl", record)
    return record


def _assemble_error_record(
    *,
    item: Dict[str, Any],
    inst: Dict[str, Any],
    adapter: Adapter,
    session_id: str,
    manifest: Dict[str, Any],
    manifest_sha: str,
    env_id: str,
    taskset_cpus: Optional[str],
    nice_val: Optional[int],
    lock_file: Optional[str],
    lock_wait_s: float,
    load_gate_info: Dict[str, Any],
    kill_s: float,
    run_dir: Path,
    session_dir: Path,
    err_msg: str,
    timestamp: str,
) -> Dict[str, Any]:
    key = item["key"]
    return {
        "schema": "dm4ri-bench-run/1",
        "run_id": item["run_id"],
        "key": key,
        "session": session_id,
        "harness": {"version": HARNESS_VERSION},
        "manifest": {"id": manifest["id"], "version": manifest["version"], "sha256": manifest_sha},
        "instance": {
            "id": inst["id"],
            "kind": inst["kind"],
            "tier": inst["tier"],
            "d_ref": inst.get("d_ref"),
            "bracket_ref": inst.get("bracket_ref"),
            "compare3_id": inst.get("compare3_id"),
        },
        "engine": adapter.describe(),
        "options": {
            "league": key["league"],
            "threads": int(item["threads"]),
            "budget_s": float(item["budget_s"]),
            "seed": int(item["seed"]),
            "argv_extra": [],
            "hints": {},
            "canary": key.get("canary"),
            "kill_after_s": round(kill_s, 3),
            "cpus": taskset_cpus,
        },
        "placement": {
            "taskset": taskset_cpus,
            "nice": nice_val,
            "exclusive": False,
            "lock": lock_file,
            "lock_wait_sec": lock_wait_s,
            "load_gate": load_gate_info,
        },
        "result": {
            "status": "error",
            "status_engine": "error",
            "lb": None,
            "ub": None,
            "exact": False,
            "wall_sec": 0.0,
            "cap_sec": float(item["budget_s"]),
            "prep_s": 0.0,
            "cpu_user_s": 0.0,
            "cpu_sys_s": 0.0,
            "max_rss_kb": 0,
            "rc": None,
            "signal": None,
            "killed": False,
            "t_internal": {"t_stop": None, "t_input": None, "t_import": None, "t_call": None},
            "stop_reason": None,
            "stop_code": None,
            "t_ub": None,
            "t_lb": None,
            "rw_steps": None,
            "size": {"full": None, "searched": None},
            "cc_rounds": [],
            "threads_used": {"requested": int(item["threads"]), "rw": None, "cc": None},
            "rw_opt": None,
            "k0": None,
            "witness": None,
            "events": [],
            "parse": "error",
            "engine_info": {},
            "error": {"kind": "harness", "msg": err_msg[:400]},
        },
        "verify": {
            "witness_weight": None,
            "witness_verified_numpy": None,
            "witness_check_detail": None,
            "ub_verified": None,
            "lb_audit": "n/a",
            "matches_ref": None,
            "violations": [],
            "unverified_reason": "harness error before execution",
        },
        "env": {
            "env_id": env_id,
            "busy_cpu_s": 0.0,
            "child_cpu_s": 0.0,
            "harness_cpu_s": 0.0,
            "ext_cpu_s": 0.0,
            "steal_s": 0.0,
            "contaminated": False,
            "started_at": timestamp,
            "ended_at": timestamp,
        },
        "raw": {
            "dir": str(run_dir.relative_to(session_dir)),
            "stdout": str((run_dir / "out.txt").relative_to(session_dir)),
            "stderr": str((run_dir / "err.txt").relative_to(session_dir)),
            "argv": [],
        },
    }


def _parse_int_csv(val: Optional[str]) -> Optional[List[int]]:
    if val is None:
        return None
    return [int(x.strip()) for x in val.split(",") if x.strip()]


def _parse_str_csv(val: Optional[str]) -> Optional[List[str]]:
    if val is None:
        return None
    return [x.strip() for x in val.split(",") if x.strip()]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="dist-m4ri engine-ready certification benchmark harness.")
    p.add_argument("--manifest", type=Path, default=MANIFEST_DIR / "certify_v1.json")
    p.add_argument("--session", type=Path, required=True, help="Session output directory (restartable).")
    p.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    p.add_argument("--cache-dir", type=Path, default=None, help="Content-addressed DEM cache directory.")
    p.add_argument(
        "--binary", type=Path, default=None, help="Path to dist_m4ri binary (default: <repo>/src/dist_m4ri)."
    )
    p.add_argument("--qec-cert-root", type=Path, default=None, help="Path to read-only qec-cert repo.")
    p.add_argument(
        "--engines",
        type=str,
        default="dist_m4ri:default",
        help="Comma-separated engine specs (e.g. dist_m4ri:default,dist_m4ri.py:off,dist_m4ri.py:auto).",
    )
    p.add_argument("--track", choices=("certify", "scaling"), default="certify")
    p.add_argument("--tiers", type=str, default=None, help="Comma-separated tiers (e.g. T0,T1).")
    p.add_argument("--instances", type=str, default=None, help="Comma-separated instance ids.")
    p.add_argument("--only-small", action="store_true", help="Skip large DEMs (recipe.large=true).")
    p.add_argument("--threads", type=str, default=None, help="Override thread list (e.g. 2 or 128,64).")
    p.add_argument("--seeds", type=str, default=None, help="Override seed list (e.g. 1,2,3; all >= 1).")
    p.add_argument("--budget", type=float, default=None, help="Override time budget in seconds.")
    p.add_argument("--order-seed", type=int, default=20261124, help="Fixed RNG seed for per-tier run permutation.")
    p.add_argument("--canary-every", type=float, default=1800.0, help="Canary interval in seconds (0 to disable).")
    p.add_argument("--canary-instance", type=str, default=None, help="Override canary instance id.")
    p.add_argument("--lock-file", type=str, default=DEFAULT_LOCK_FILE, help="flock path (shared with qec-cert).")
    p.add_argument("--no-lock", action="store_true", help="Test flag: do not acquire flock (sets exclusive=false).")
    p.add_argument("--max-busy-cores", type=float, default=2.0, help="Load gate threshold in busy cores.")
    p.add_argument("--load-window", type=float, default=2.0, help="Load gate /proc/stat sampling window (s).")
    p.add_argument("--load-poll", type=float, default=30.0, help="Seconds between load gate retries.")
    p.add_argument("--load-max-wait", type=float, default=None, help="Optional max seconds to wait at load gate.")
    p.add_argument("--no-load-gate", action="store_true", help="Test flag: skip load gate (for pinned T0 tests).")
    p.add_argument("--taskset", type=str, default=None, help="Override CPU affinity list (e.g. 126-127).")
    p.add_argument("--nice", type=int, default=None, help="Optional nice increment (e.g. 19 for tests).")
    p.add_argument("--retry", action="store_true", help="Retry prior runs with error.kind == 'harness'.")
    p.add_argument("--retry-all-errors", action="store_true", help="Retry all prior error runs.")
    p.add_argument("--plan-only", "--dry-run", action="store_true", dest="plan_only", help="Write plan and exit.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    repo_root = args.repo_root.resolve()
    cache_dir = (args.cache_dir if args.cache_dir is not None else mi.default_cache_dir(repo_root)).resolve()
    binary_path = (args.binary if args.binary is not None else (repo_root / "src" / "dist_m4ri")).resolve()
    session_dir = args.session.resolve()
    session_dir.mkdir(parents=True, exist_ok=True)
    log_path = session_dir / "harness.log"

    def log(msg: str) -> None:
        line = f"[{utc_now_iso()}] {msg}"
        print(line, flush=True)
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(line + "\n")

    manifest, manifest_sha = mi.load_manifest(args.manifest.resolve())
    by_id = {inst["id"]: inst for inst in manifest["instances"]}

    engine_specs = _parse_str_csv(args.engines) or ["dist_m4ri:default"]
    adapters: List[Adapter] = [
        build_adapter(
            spec,
            repo_root=repo_root,
            binary_path=binary_path,
            cache_dir=cache_dir,
            qec_cert_root=args.qec_cert_root,
        )
        for spec in engine_specs
    ]
    adapter_map: Dict[Tuple[str, str], Adapter] = {(a.name, a.config): a for a in adapters}

    plan = build_plan(
        manifest,
        adapters,
        track=args.track,
        tiers=_parse_str_csv(args.tiers),
        instance_ids=_parse_str_csv(args.instances),
        only_small=args.only_small,
        threads_override=_parse_int_csv(args.threads),
        seeds_override=_parse_int_csv(args.seeds),
        budget_override=args.budget,
        order_seed=args.order_seed,
    )
    atomic_write_json(session_dir / "plan.json", plan)
    log(
        f"Plan written: track={args.track}, runs={plan['total_runs']}, "
        f"engines={plan['options']['engines']}, order_seed={args.order_seed}"
    )

    # Verify all required inputs (plus canary instance if enabled) and warm page cache
    needed_ids = {r["instance"] for r in plan["runs"]}
    canary_id = args.canary_instance or manifest.get("canary_instance")
    if args.canary_every > 0 and canary_id and canary_id in by_id and not args.plan_only:
        needed_ids.add(canary_id)
    needed_insts = [by_id[iid] for iid in sorted(needed_ids) if iid in by_id]
    inputs_summary = mi.verify_manifest_inputs(
        manifest, needed_insts, repo_root=repo_root, cache_dir=cache_dir, warm_page_cache=not args.plan_only
    )
    log(
        f"Inputs verified: {inputs_summary['files']} files "
        f"({inputs_summary['bytes'] / 1e6:.2f} MB), sha256_ok=True, "
        f"page_cache_warmed={inputs_summary['page_cache_warmed']}"
    )

    env_obj, env_id = collect_env_json(
        repo_root=repo_root,
        binary_path=binary_path,
        manifest=manifest,
        manifest_sha=manifest_sha,
        inputs_summary=inputs_summary,
        qec_cert_root=args.qec_cert_root,
    )
    env_path = session_dir / "env.json"
    if env_path.is_file():
        prev_env = json.loads(env_path.read_text(encoding="utf-8"))
        prev_cmp = {k: v for k, v in prev_env.items() if k != "created"}
        curr_cmp = {k: v for k, v in env_obj.items() if k != "created"}
        if prev_cmp == curr_cmp:
            env_obj = prev_env
            env_id = hashlib.sha1(canonical_json(env_obj).encode("utf-8")).hexdigest()
        else:
            prev_id = hashlib.sha1(canonical_json(prev_env).encode("utf-8")).hexdigest()
            atomic_write_json(session_dir / "envs" / f"{prev_id}.json", prev_env)
            atomic_write_json(env_path, env_obj)
            atomic_write_json(session_dir / "envs" / f"{env_id}.json", env_obj)
    else:
        atomic_write_json(env_path, env_obj)

    session_meta_path = session_dir / "session.json"
    if session_meta_path.is_file():
        session_meta = json.loads(session_meta_path.read_text(encoding="utf-8"))
        session_id = session_meta["session"]
    else:
        session_id = f"{utc_now_iso()}-{socket.gethostname()}"
        session_meta = {
            "schema": "dm4ri-bench-session/1",
            "session": session_id,
            "created": utc_now_iso(),
            "harness_version": HARNESS_VERSION,
            "manifest": {"id": manifest["id"], "version": manifest["version"], "sha256": manifest_sha},
            "env_id": env_id,
            "argv": sys.argv,
        }
        atomic_write_json(session_meta_path, session_meta)

    if args.plan_only:
        for idx, item in enumerate(plan["runs"][:30], 1):
            print(
                f"  {idx:3d}. [{item['tier']}] {item['instance']} "
                f"{item['engine']}:{item['config']}@p{item['threads']} "
                f"T={item['budget_s']:g}s s={item['seed']} ({item['run_id'][:8]})"
            )
        if len(plan["runs"]) > 30:
            print(f"  ... ({len(plan['runs']) - 30} more runs in {session_dir / 'plan.json'})")
        return 0

    results_path = session_dir / "results.jsonl"
    existing = load_existing_records(results_path)
    verifier = vf.VerifierContext(repo_root=repo_root, cache_dir=cache_dir)
    lock_file = None if args.no_lock else args.lock_file
    use_load_gate = not args.no_load_gate

    canary_adapter: Optional[Adapter] = None
    canary_count = sum(
        1 for rec in existing.values() if (rec.get("options") or {}).get("canary") is not None
    )
    last_canary_t = 0.0
    if args.canary_every > 0 and canary_id and canary_id in by_id:
        canary_adapter = adapter_map.get(("dist_m4ri", "default")) or build_adapter(
            "dist_m4ri:default", repo_root=repo_root, binary_path=binary_path, cache_dir=cache_dir
        )

    def maybe_run_canary() -> None:
        nonlocal canary_count, last_canary_t
        if canary_adapter is None or not canary_id or canary_id not in by_id:
            return
        now_m = time.monotonic()
        if last_canary_t > 0.0 and (now_m - last_canary_t) < args.canary_every:
            return
        c_inst = by_id[canary_id]
        c_tier = manifest["tiers"][c_inst["tier"]]
        c_thr = int((_parse_int_csv(args.threads) or [c_tier["threads"][0]])[0])
        c_budget = float(args.budget) if args.budget is not None else float(c_tier["budget_s"])
        canary_count += 1
        c_key = {
            "manifest": f"{manifest['id']}/{manifest['version']}",
            "instance": c_inst["id"],
            "engine": canary_adapter.name,
            "config": canary_adapter.config,
            "league": "default",
            "threads": c_thr,
            "budget_s": c_budget,
            "seed": 1,
            "canary": canary_count,
            "harness_major": HARNESS_MAJOR,
        }
        c_item = {
            "run_id": compute_run_id(c_key),
            "key": c_key,
            "tier": c_inst["tier"],
            "instance": c_inst["id"],
            "engine": canary_adapter.name,
            "config": canary_adapter.config,
            "league": "default",
            "threads": c_thr,
            "budget_s": c_budget,
            "seed": 1,
        }
        rec = execute_single_run(
            item=c_item,
            inst=c_inst,
            adapter=canary_adapter,
            verifier=verifier,
            session_dir=session_dir,
            session_id=session_id,
            manifest=manifest,
            manifest_sha=manifest_sha,
            env_id=env_id,
            lock_file=lock_file,
            use_load_gate=use_load_gate,
            max_busy_cores=args.max_busy_cores,
            load_window_s=args.load_window,
            load_poll_s=args.load_poll,
            load_max_wait_s=args.load_max_wait,
            taskset_override=args.taskset,
            nice_val=args.nice,
            log_fn=log,
        )
        existing[c_item["run_id"]] = rec
        last_canary_t = time.monotonic()
        res = rec["result"]
        log(
            f"CANARY #{canary_count} {c_inst['id']} p={c_thr}: status={res['status']} "
            f"[{res['lb']}, {res['ub']}] wall={res['wall_sec']:.3f}s ext_cpu={rec['env']['ext_cpu_s']:.2f}s"
        )

    skipped = 0
    executed = 0
    wrong_cnt = 0
    for idx, item in enumerate(plan["runs"], 1):
        rid = item["run_id"]
        if should_skip_run(existing.get(rid), retry=args.retry, retry_all_errors=args.retry_all_errors):
            skipped += 1
            continue
        maybe_run_canary()
        inst = by_id[item["instance"]]
        ad = adapter_map[(item["engine"], item["config"])]
        rec = execute_single_run(
            item=item,
            inst=inst,
            adapter=ad,
            verifier=verifier,
            session_dir=session_dir,
            session_id=session_id,
            manifest=manifest,
            manifest_sha=manifest_sha,
            env_id=env_id,
            lock_file=lock_file,
            use_load_gate=use_load_gate,
            max_busy_cores=args.max_busy_cores,
            load_window_s=args.load_window,
            load_poll_s=args.load_poll,
            load_max_wait_s=args.load_max_wait,
            taskset_override=args.taskset,
            nice_val=args.nice,
            log_fn=log,
        )
        existing[rid] = rec
        executed += 1
        res = rec["result"]
        ver = rec["verify"]
        if res["status"] == "wrong":
            wrong_cnt += 1
        log(
            f"[{idx}/{plan['total_runs']}] {item['instance']} {item['engine']}:{item['config']}@p{item['threads']} "
            f"s={item['seed']}: {res['status']} [{res['lb']}, {res['ub']}] "
            f"wall={res['wall_sec']:.3f}s rss={res['max_rss_kb']}kB "
            f"wit_ok={ver.get('witness_verified_numpy')} ext={rec['env']['ext_cpu_s']:.2f}s"
        )

    audit = vf.audit_records(list(existing.values()), manifest)
    atomic_write_json(session_dir / "audit.json", audit)
    log(
        f"Session finished: executed={executed}, skipped={skipped}, total={len(existing)}, "
        f"wrong={len(audit['wrong_runs'])}, lb_conflicts={len(audit['lb_conflicts'])}, "
        f"contaminated={len(audit['contaminated_runs'])}"
    )
    return 0 if audit["ok"] and wrong_cnt == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
