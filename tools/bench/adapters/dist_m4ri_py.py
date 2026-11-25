"""Adapter for the Python `dist_m4ri.py` library (`compute_*_distance`).

Configurations:
  - config='off'  (engine=None):   Python wrapper over the `dist_m4ri` binary (cache disabled).
  - config='auto' (engine='auto'): Two-stage portfolio (`dist_m4ri` -> `qec_cert`).
  - config='cert' (engine='cert'): Direct delegation to `qec_cert` through `dist_m4ri.py`.

Each run executes in a fresh Python subprocess (`py_worker.py`) so process wall time includes
imports, Python wrapper overhead, and (for circuits) Stim -> DEM conversion.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
MANIFEST_DIR = REPO_ROOT / "benchmark" / "manifest"
if str(MANIFEST_DIR) not in sys.path:
    sys.path.insert(0, str(MANIFEST_DIR))
import make_inputs as mi  # noqa: E402

from . import Adapter, git_info  # noqa: E402
from .dist_m4ri_bin import DistM4riBinAdapter, parse_binary_stderr  # noqa: E402


WORKER_SCRIPT = Path(__file__).resolve().parent / "py_worker.py"


def find_qec_cert_root(repo_root: Path, explicit: Optional[Path] = None) -> Optional[Path]:
    """Locate the read-only qec-cert repository if present."""
    candidates = []
    if explicit is not None:
        candidates.append(explicit)
    env_repo = os.environ.get("QEC_CERT_REPO")
    if env_repo:
        candidates.append(Path(env_repo))
    candidates.append(repo_root.parent / "qec-cert")
    try:
        spec = importlib.util.find_spec("qec_cert")
        if spec and spec.origin:
            candidates.append(Path(spec.origin).resolve().parent.parent)
    except Exception:
        pass
    for cand in candidates:
        if cand and (cand / "qec_cert").is_dir():
            return cand.resolve()
    return None


def find_qec_cert_lib(qec_root: Optional[Path]) -> Optional[Path]:
    """Locate libqeccert.so if present."""
    env_lib = os.environ.get("QEC_CERT_LIB")
    if env_lib and Path(env_lib).is_file():
        return Path(env_lib).resolve()
    if qec_root is not None:
        for rel in ("qec_cert/libqeccert.so", "build/libqeccert.so"):
            cand = qec_root / rel
            if cand.is_file():
                return cand.resolve()
    return None


def _py_version_of(mod_path: Path) -> Optional[str]:
    if not mod_path.is_file():
        return None
    for line in mod_path.read_text(encoding="utf-8").splitlines()[:80]:
        m = re.match(r'^__version__\s*=\s*["\']([^"\']+)["\']', line.strip())
        if m:
            return m.group(1)
    return None


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


class DistM4riPyAdapter(Adapter):
    """Adapter for `dist_m4ri.py` with `engine` in {None/'off', 'auto', 'cert'}."""

    name = "dist_m4ri.py"
    interface = "python"

    def __init__(
        self,
        config: Optional[str] = "off",
        repo_root: Path = REPO_ROOT,
        binary_path: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        python_exe: Optional[str] = None,
        qec_cert_root: Optional[Path] = None,
    ) -> None:
        norm = "off" if config in (None, "", "none", "off", "default") else str(config)
        if norm not in ("off", "auto", "cert"):
            raise ValueError(f"Unsupported dist_m4ri.py config {config!r} (expected 'off', 'auto', or 'cert')")
        self.config = norm
        self.engine_arg: Optional[str] = None if norm == "off" else norm
        self.repo_root = repo_root.resolve()
        self.py_path = self.repo_root / "dist_m4ri.py"
        self.binary = (binary_path if binary_path is not None else (self.repo_root / "src" / "dist_m4ri")).resolve()
        self.cache_dir = (cache_dir if cache_dir is not None else mi.default_cache_dir(self.repo_root)).resolve()
        self.python_exe = python_exe or sys.executable
        self.qec_cert_root = find_qec_cert_root(self.repo_root, explicit=qec_cert_root)
        self._bin_adapter = DistM4riBinAdapter(
            binary_path=self.binary, config="default", repo_root=self.repo_root, cache_dir=self.cache_dir
        )
        self._desc_cache: Optional[Dict[str, Any]] = None
        self._last_threads: Optional[int] = None
        self._last_budget_s: Optional[float] = None

    def describe(self) -> Dict[str, Any]:
        if self._desc_cache is not None:
            return dict(self._desc_cache)
        commit, dirty = git_info(self.repo_root)
        bin_desc = self._bin_adapter.describe()
        deps: Dict[str, Any] = {
            "python": sys.version.split()[0],
            "stim": _pkg_version("stim"),
            "numpy": _pkg_version("numpy"),
            "scipy": _pkg_version("scipy"),
            "dist_m4ri_bin_version": bin_desc.get("version"),
            "dist_m4ri_bin_sha256": bin_desc.get("sha256"),
        }
        if self.qec_cert_root is not None:
            q_commit, q_dirty = git_info(self.qec_cert_root)
            q_lib = find_qec_cert_lib(self.qec_cert_root)
            deps["qec_cert"] = {
                "commit": q_commit,
                "dirty": q_dirty,
                "lib_sha256": mi.sha256_file(q_lib) if q_lib else None,
            }
        self._desc_cache = {
            "name": self.name,
            "interface": self.interface,
            "config": self.config,
            "version": _py_version_of(self.py_path),
            "commit": commit,
            "dirty": dirty,
            "sha256": mi.sha256_file(self.py_path) if self.py_path.is_file() else None,
            "deps": deps,
            "features": {"engine": self.engine_arg, "cache_disabled": True},
        }
        return dict(self._desc_cache)

    def supports(self, inst: Dict[str, Any]) -> bool:
        if not self.py_path.is_file():
            return False
        if inst["kind"] not in ("css_sector", "classical", "hl", "dem", "circuit"):
            return False
        if self.config == "cert":
            return find_qec_cert_lib(self.qec_cert_root) is not None
        return self.binary.is_file()

    def prepare(self, inst: Dict[str, Any]) -> Dict[str, Any]:
        files = mi.resolve_instance_files(
            inst, repo_root=self.repo_root, cache_dir=self.cache_dir, verify_sha=False
        )
        dem_mode = "stim_auto" if inst["kind"] == "circuit" else "canonical"
        return {
            "prep_s": 0.0,
            "files": {k: str(v) for k, v in files.items()},
            "dem_mode": dem_mode,
        }

    def command(
        self,
        inst: Dict[str, Any],
        threads: int,
        budget_s: float,
        seed: int,
        cfg: Dict[str, Any],
    ) -> Dict[str, Any]:
        if seed < 1:
            raise ValueError(f"Harness rule m9: seed must be >= 1 (got {seed})")
        prep = cfg.get("prep") or self.prepare(inst)
        d_ref = inst.get("d_ref")
        br = inst.get("bracket_ref")
        may_exceed_25 = (
            (d_ref is not None and int(d_ref) > 25)
            or (br is not None and int(br[1]) > 25)
            or (d_ref is None and br is None)
        )
        # Forward DBG_SUMMARY|DBG_PROGRESS|DBG_CODEWORDS (and DBG_LARGE if d>25) when binary may run
        dbg = 0 if self.config == "cert" else ((35 | 2048) if may_exceed_25 else 35)
        spec = {
            "repo_root": str(self.repo_root),
            "binary_path": str(self.binary),
            "qec_cert_root": str(self.qec_cert_root) if self.qec_cert_root else None,
            "kind": inst["kind"],
            "files": prep["files"],
            "engine": self.engine_arg,
            "threads": int(threads),
            "budget_s": float(budget_s),
            "seed": int(seed),
            "debug": dbg,
        }
        self._last_threads = threads
        self._last_budget_s = float(budget_s)
        argv = [self.python_exe, str(WORKER_SCRIPT), json.dumps(spec)]
        argv_extra = [
            f"engine={self.engine_arg!r}",
            f"threads={threads}",
            f"timeout={budget_s:g}",
            f"seed={seed}",
            f"debug={dbg}",
        ]
        return {"argv": argv, "argv_extra": argv_extra, "env": {}, "cwd": str(self.repo_root)}

    def parse(
        self,
        rc: Optional[int],
        stdout: str,
        stderr: str,
        *,
        killed: bool = False,
        prep: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        dem_mode = (prep or {}).get("dem_mode", "canonical")
        py_obj: Optional[Dict[str, Any]] = None
        for line in reversed((stdout or "").splitlines()):
            if line.startswith("@@PYRESULT "):
                try:
                    py_obj = json.loads(line[len("@@PYRESULT "):].strip())
                    break
                except ValueError:
                    pass

        # Parse any underlying binary progress/witness output emitted to stdout or stderr
        bin_parsed = parse_binary_stderr(
            "",
            (stdout or "") + "\n" + (stderr or ""),
            rc=0,
            killed=killed,
            config="default",
            dem_mode=dem_mode,
            requested_threads=self._last_threads,
        )

        if killed:
            return {
                "status": "timeout",
                "status_engine": "timeout",
                "lb": bin_parsed.get("lb"),
                "ub": bin_parsed.get("ub"),
                "exact": False,
                "t_internal": {
                    "t_stop": bin_parsed["t_internal"].get("t_stop"),
                    "t_input": None,
                    "t_import": None,
                    "t_call": None,
                },
                "stop_reason": "killed by harness at kill_after_s",
                "stop_code": "killed",
                "t_ub": bin_parsed.get("t_ub"),
                "t_lb": bin_parsed.get("t_lb"),
                "rw_steps": bin_parsed.get("rw_steps"),
                "size": bin_parsed.get("size", {"full": None, "searched": None}),
                "cc_rounds": bin_parsed.get("cc_rounds", []),
                "threads_used": bin_parsed.get(
                    "threads_used", {"requested": self._last_threads, "rw": None, "cc": None}
                ),
                "rw_opt": bin_parsed.get("rw_opt"),
                "k0": bin_parsed.get("k0"),
                "witness": bin_parsed.get("witness"),
                "events": bin_parsed.get("events", []),
                "parse": "worker",
                "engine_info": {},
                "error": None,
            }

        if py_obj is None:
            msg = (stderr or stdout or "").strip().splitlines()
            return {
                "status": "error",
                "status_engine": "error",
                "lb": None,
                "ub": None,
                "exact": False,
                "t_internal": {"t_stop": None, "t_input": None, "t_import": None, "t_call": None},
                "stop_reason": None,
                "stop_code": None,
                "t_ub": None,
                "t_lb": None,
                "rw_steps": None,
                "size": {"full": None, "searched": None},
                "cc_rounds": [],
                "threads_used": {"requested": self._last_threads, "rw": None, "cc": None},
                "rw_opt": None,
                "k0": None,
                "witness": None,
                "events": [],
                "parse": "worker",
                "engine_info": {},
                "error": {"kind": "engine", "msg": (msg[-1][:240] if msg else f"exit code {rc}")},
            }

        if py_obj.get("error"):
            return {
                "status": "error",
                "status_engine": "error",
                "lb": None,
                "ub": None,
                "exact": False,
                "t_internal": {
                    "t_stop": None,
                    "t_input": None,
                    "t_import": py_obj.get("t_import"),
                    "t_call": py_obj.get("t_call"),
                },
                "stop_reason": None,
                "stop_code": None,
                "t_ub": None,
                "t_lb": None,
                "rw_steps": None,
                "size": {"full": None, "searched": None},
                "cc_rounds": [],
                "threads_used": {"requested": self._last_threads, "rw": None, "cc": None},
                "rw_opt": None,
                "k0": None,
                "witness": None,
                "events": [],
                "parse": "worker",
                "engine_info": {"qec_imported": py_obj.get("qec_imported")},
                "error": {"kind": "engine", "msg": str(py_obj["error"])},
            }

        d_info = py_obj.get("d_info") or [0, 0, 0]
        lb = int(d_info[0])
        ub: Optional[int] = int(d_info[1]) if int(d_info[1]) > 0 else None
        rw_steps = int(d_info[2])
        stats = py_obj.get("last_run_stats") or {}
        k0 = bool(stats.get("k0")) or (lb == 0 and int(d_info[1]) == 0 and py_obj.get("dist") == 0)
        if k0:
            ub = 0

        exact = bool((lb > 0 and ub is not None and lb == ub) or k0)
        auto_info = stats.get("engine_auto") or {}
        cert_info = stats.get("engine_cert") or {}
        stage_b_ran = auto_info.get("stage_b") is not None
        stage_c_ran = auto_info.get("stage_c") is not None
        ended_in_qec = (self.config == "cert") or (stage_b_ran and not stage_c_ran)

        budget_s = py_obj.get("budget_s")
        if budget_s is None:
            budget_s = self._last_budget_s
        budget_f = float(budget_s) if budget_s is not None else None
        t_call = py_obj.get("t_call")
        t_call_f = float(t_call) if t_call is not None else None
        cert_t = cert_info.get("time")
        cert_t_f = float(cert_t) if cert_t is not None else None
        auto_reason = str(auto_info.get("reason") or "")
        auto_reason_lo = auto_reason.lower()

        # M1: Do not inherit Stage A's sub-budget stop_code/stop_reason/t_stop when Stage B or cert finished.
        if ended_in_qec:
            stop_code: Optional[str] = None
            stop_reason: Optional[str] = None
            t_stop_val: Optional[float] = None
            t_input_val: Optional[float] = None
        else:
            stop_code = bin_parsed.get("stop_code")
            stop_reason = bin_parsed.get("stop_reason")
            t_stop_val = bin_parsed["t_internal"].get("t_stop")
            t_input_val = bin_parsed["t_internal"].get("t_input")

        if exact:
            status = "certified"
            if k0:
                stop_code = "k0"
            elif ended_in_qec or stop_code in (None, "timeout"):
                stop_code = "bounds"
            if ended_in_qec:
                default_sr = "certified by qec_cert" if self.config == "cert" else "exact after Stage B"
                stop_reason = auto_reason or default_sr
        else:
            # Robust timeout detection for cert, auto Stage B (incl. "0 s left after qec_cert/Stage B"), and binary.
            near_budget = bool(
                budget_f is not None
                and budget_f > 0.0
                and (
                    (t_call_f is not None and t_call_f >= 0.95 * budget_f)
                    or (cert_t_f is not None and cert_t_f >= 0.95 * budget_f)
                )
            )
            stage_b_exhausted = bool(
                "s left after" in auto_reason_lo
                or "timeout" in auto_reason_lo
                or (
                    budget_f is not None
                    and budget_f > 0.0
                    and auto_info.get("stage_b_time") is not None
                    and (
                        float(auto_info.get("stage_a_time") or 0.0)
                        + float(auto_info.get("stage_b_time") or 0.0)
                    )
                    >= 0.95 * budget_f
                )
            )
            if stop_code == "timeout" or stage_b_exhausted or near_budget:
                status = "timeout"
                stop_code = "timeout"
                if not stop_reason:
                    stop_reason = auto_reason or (
                        f"budget exhausted ({budget_f:g}s)" if budget_f is not None else "timeout"
                    )
            else:
                status = "bracket"
                if ended_in_qec and not stop_reason:
                    stop_reason = auto_reason or "qec_cert ended with open bracket"

        # M2: Attach binary witness ONLY when its weight equals the final ub.
        raw_wit = bin_parsed.get("witness")
        if raw_wit is not None and ub is not None and int(raw_wit.get("weight", -1)) == int(ub):
            witness: Optional[Dict[str, Any]] = raw_wit
        else:
            witness = None

        engine_info: Dict[str, Any] = {
            "qec_imported": py_obj.get("qec_imported"),
            "engine_auto": stats.get("engine_auto"),
            "engine_cert": stats.get("engine_cert"),
        }
        if witness is None and ub is not None and ub > 0:
            if self.config == "cert" or stage_b_ran:
                engine_info["witness_unavailable_reason"] = (
                    "dist_m4ri.py engine='cert'/'auto' Stage B does not expose qec_cert witness support"
                )

        return {
            "status": status,
            "status_engine": status,
            "lb": lb,
            "ub": ub,
            "exact": exact,
            "t_internal": {
                "t_stop": t_stop_val,
                "t_input": t_input_val,
                "t_import": py_obj.get("t_import"),
                "t_call": py_obj.get("t_call"),
            },
            "stop_reason": stop_reason,
            "stop_code": stop_code,
            "t_ub": bin_parsed.get("t_ub"),
            "t_lb": bin_parsed.get("t_lb"),
            "rw_steps": rw_steps,
            "size": bin_parsed.get("size", {"full": None, "searched": None}),
            "cc_rounds": bin_parsed.get("cc_rounds", []),
            "threads_used": bin_parsed.get("threads_used", {"requested": self._last_threads, "rw": None, "cc": None}),
            "rw_opt": bin_parsed.get("rw_opt"),
            "k0": k0,
            "witness": witness,
            "events": bin_parsed.get("events", []),
            "parse": "worker",
            "engine_info": engine_info,
            "error": None,
        }
