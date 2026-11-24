"""Adapter for the native `dist_m4ri` C binary.

Supports:
  - config='default': method=3 with ONLY threads, timeout, and seed (plus opt-in json=1 stamp=1 debug=3
    when the binary supports it, or debug=35 / debug=2083 for full witness extraction on older binaries).
  - config='cc_fixed': Track S thread scaling (method=2 wmax=d_ref-1 fixed work).

Auto-detects `@binstamp` machine-readable output by probing `<binary> --help json` once per binary.
Falls back cleanly to stdout + stderr parsing (with line-interleaving recovery) for older binaries
and for runs killed by the harness before `# RESULT` is emitted.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
MANIFEST_DIR = REPO_ROOT / "benchmark" / "manifest"
if str(MANIFEST_DIR) not in sys.path:
    sys.path.insert(0, str(MANIFEST_DIR))
import make_inputs as mi  # noqa: E402

from . import Adapter, git_info  # noqa: E402


RE_RESULT_LINE = re.compile(r"^\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$")
RE_STAMP = re.compile(r"\s+\[t=(\d+(?:\.\d+)?)s\]\s*$")
RE_STOP = re.compile(r"^#\s*stopped after ([0-9.]+)s:\s*(.*)$")
RE_STOP_NOSEARCH = re.compile(r"^#\s*stopped:\s*(.*)$")
RE_UB = re.compile(r"RW(?:\s*\(ksub=\d+\))?\s+found new upper bound cw of weight (\d+)")
RE_CC_ROUND = re.compile(
    r"^#\s*CC round w=(\d+) (?:finished|completed) in ([0-9.]+)s"
    r"(?:\s*\((\d+) CC threads,\s*(\d+) RW threads\))?"
)
RE_DMIN_ARROW = re.compile(r"->\s*dmin=(\d+)")
RE_CC_WORK = re.compile(
    r"^#\s*CC round w=(\d+) work:\s*([0-9.eE+-]+) thread-s measured.*?(\d+) nodes"
)
RE_CW = re.compile(
    r"^#\s*(?:\[thread \d+\]\s*)?(RW|CC|finC:?)\s+(?:lightest\s+)?codeword of weight (\d+) "
    r"\(1-based columns\):\s*(.*)$"
)
RE_H_SIZE = re.compile(r"^#\s+H:\s*(\d+)\s*x\s*(\d+)")
RE_MERGED = re.compile(r"merged \d+ equal columns of the DEM:\s*(\d+) of (\d+) columns kept")
RE_DROPPED = re.compile(r"dropped \d+ of (\d+) detectors")
RE_RW_INFO = re.compile(
    r"RW information sets:\s*n=(\d+),\s*rank\(H\)=(\d+),\s*steps=(\d+).*?(\d+) RW threads"
)
RE_RW_OPT = re.compile(r"rw_opt=(0x[0-9a-fA-F]+)")
RE_RUN_PLAN_THREADS = re.compile(r"^#\s*run plan:.*?,\s*(\d+) threads \(RW on up to (\d+)\)")
RE_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def _strip_ansi(text: str) -> str:
    return RE_ANSI.sub("", text or "")


def _unscramble_stderr(stderr: str) -> List[str]:
    """Recover stderr lines when two threads interleave `# ...` messages mid-line."""
    raw_lines = stderr.splitlines()
    out: List[str] = []
    pending: Optional[str] = None
    for line in raw_lines:
        if pending is not None and not line.startswith("#"):
            out.append(pending + line)
            pending = None
            continue
        if pending is not None:
            out.append(pending)
            pending = None
        if line.startswith("#"):
            # Look for an embedded '# ' starting a concurrent log line after position 2
            idx = line.find("# ", 2)
            if idx > 0 and ("codeword of weight" in line[:idx] or "found new upper bound" in line[:idx]):
                pending = line[:idx]
                out.append(line[idx:])
                continue
        out.append(line)
    if pending is not None:
        out.append(pending)
    return out


def _map_stop_code(text: Optional[str]) -> Optional[str]:
    if not text:
        return None
    low = text.lower()
    if "bounds coincide" in low:
        return "bounds"
    if "cc found" in low:
        return "cc_found"
    if "timeout" in low:
        return "timeout"
    if "rw convergence" in low:
        return "rw_converged"
    if "all rw steps" in low or "rw steps done" in low:
        return "rw_steps"
    if "predicted" in low or "cc_slow" in low or "not to finish before" in low:
        return "cc_slow"
    if "wmax" in low or "largest weight" in low or "cc done" in low:
        return "cc_done"
    if "wmin" in low:
        return "wmin"
    if "maxc" in low:
        return "maxc"
    if "collected" in low:
        return "collected"
    if "dstop" in low:
        return "dstop"
    if "k = 0" in low or "k=0" in low or "no logical" in low:
        return "k0"
    return "other"


def _parse_stamp(line: str) -> Tuple[str, Optional[float]]:
    m = RE_STAMP.search(line)
    if not m:
        return line, None
    return line[:m.start()], float(m.group(1))


def parse_binary_stderr(
    stdout: str,
    stderr: str,
    *,
    rc: Optional[int] = 0,
    killed: bool = False,
    config: str = "default",
    dem_mode: str = "canonical",
    requested_threads: Optional[int] = None,
) -> Dict[str, Any]:
    """Parse stdout and stderr of `dist_m4ri` (using `# RESULT {json}` if present, else stderr text)."""
    clean_err = _strip_ansi(stderr or "")
    lines = _unscramble_stderr(clean_err)

    # 1. Check for machine-readable `# RESULT {json}` from @binstamp
    json_obj: Optional[Dict[str, Any]] = None
    for line in reversed(lines):
        if line.startswith("# RESULT "):
            payload = line[len("# RESULT "):].strip()
            try:
                parsed = json.loads(payload)
                if isinstance(parsed, dict) and parsed.get("schema") == "dm4ri-result/1":
                    json_obj = parsed
                    break
            except ValueError:
                pass

    if json_obj is not None:
        return _from_json_result(json_obj, rc=rc, killed=killed, config=config, dem_mode=dem_mode)

    # 2. Fallback: parse stdout result line + stderr progress/summary/witness lines
    lb: Optional[int] = None
    ub: Optional[int] = None
    rw_steps: Optional[int] = None
    for line in reversed((stdout or "").splitlines()):
        m = RE_RESULT_LINE.match(line)
        if m:
            dmin_v, dmax_v, steps_v = int(m.group(1)), int(m.group(2)), int(m.group(3))
            lb = dmin_v
            ub = dmax_v if dmax_v > 0 else None
            rw_steps = steps_v
            break

    t_stop: Optional[float] = None
    stop_reason: Optional[str] = None
    stop_code: Optional[str] = None
    t_ub: Optional[float] = None
    t_lb: Optional[float] = None
    events: List[List[Any]] = []
    cc_rounds_by_w: Dict[int, Dict[str, Any]] = {}
    searched_size: Optional[List[int]] = None
    full_m: Optional[int] = None
    full_n: Optional[int] = None
    rw_opt: Optional[int] = None
    rw_threads: Optional[int] = None
    cc_threads: Optional[int] = None
    k0 = False
    prog_lb: Optional[int] = None
    prog_ub: Optional[int] = None
    witnesses: List[Dict[str, Any]] = []

    for raw_line in lines:
        line, stamp_t = _parse_stamp(raw_line)

        if "k = 0" in line or "no logical operators (k=0)" in line:
            k0 = True

        m = RE_H_SIZE.match(line)
        if m and searched_size is None:
            searched_size = [int(m.group(1)), int(m.group(2))]

        m = RE_MERGED.search(line)
        if m:
            full_n = int(m.group(2))

        m = RE_DROPPED.search(line)
        if m:
            full_m = int(m.group(1))

        m = RE_RW_OPT.search(line)
        if m:
            rw_opt = int(m.group(1), 16)

        m = RE_RUN_PLAN_THREADS.match(line)
        if m:
            cc_threads = int(m.group(1))
            if rw_threads is None:
                rw_threads = int(m.group(2))

        m = RE_RW_INFO.search(line)
        if m:
            rw_threads = int(m.group(4))
            if rw_steps is None:
                rw_steps = int(m.group(3))

        m = RE_UB.search(line)
        if m:
            w = int(m.group(1))
            if prog_ub is None or w < prog_ub:
                prog_ub = w
                if stamp_t is not None:
                    t_ub = stamp_t
                    events.append([stamp_t, prog_lb, prog_ub])

        m = RE_CC_ROUND.match(line)
        if m:
            w = int(m.group(1))
            wall_s = float(m.group(2))
            if m.group(3):
                cc_threads = max(cc_threads or 0, int(m.group(3)))
            entry = cc_rounds_by_w.setdefault(
                w, {"w": w, "thread_s": None, "wall_s": wall_s, "nodes": None, "complete": True}
            )
            entry["wall_s"] = wall_s
            entry["complete"] = True
            md = RE_DMIN_ARROW.search(line)
            if md:
                new_lb = int(md.group(1))
                if prog_lb is None or new_lb > prog_lb:
                    prog_lb = new_lb
                    if stamp_t is not None:
                        t_lb = stamp_t
                        events.append([stamp_t, prog_lb, prog_ub])

        m = RE_CC_WORK.match(line)
        if m:
            w = int(m.group(1))
            entry = cc_rounds_by_w.setdefault(
                w, {"w": w, "thread_s": None, "wall_s": None, "nodes": None, "complete": True}
            )
            entry["thread_s"] = float(m.group(2))
            entry["nodes"] = int(m.group(3))

        m = RE_CW.match(line)
        if m:
            src_raw = m.group(1).rstrip(":").lower()
            wt = int(m.group(2))
            rest = m.group(3)
            truncated = "..." in rest
            body = rest.split("...", 1)[0]
            nums = [int(tok) - 1 for tok in re.findall(r"\d+", body)]
            complete = (not truncated) and (len(nums) == wt)
            witnesses.append({
                "weight": wt,
                "support": nums,
                "base": 0,
                "source": src_raw,
                "complete": complete,
                "dem_mode": dem_mode,
            })

        m = RE_STOP.match(line)
        if m:
            t_stop = float(m.group(1))
            stop_reason = m.group(2).strip()
            stop_code = _map_stop_code(stop_reason)
        else:
            m2 = RE_STOP_NOSEARCH.match(line)
            if m2 and stop_reason is None:
                stop_reason = m2.group(1).strip()
                stop_code = _map_stop_code(stop_reason)

    # For killed runs without stdout result line, use last progress bounds
    if lb is None and prog_lb is not None:
        lb = prog_lb
    if ub is None and prog_ub is not None:
        ub = prog_ub
    if k0 and lb == 0 and ub is None:
        ub = 0

    # Pick the best witness matching ub (preferring complete support of weight == ub, else lightest)
    chosen_witness: Optional[Dict[str, Any]] = None
    if witnesses:
        if ub is not None and ub > 0:
            matching = [w for w in witnesses if w["weight"] == ub]
            if matching:
                complete_m = [w for w in matching if w["complete"]]
                chosen_witness = complete_m[-1] if complete_m else matching[-1]
        if chosen_witness is None:
            witnesses.sort(key=lambda w: (w["weight"], not w["complete"]))
            chosen_witness = witnesses[0]

    # Full size defaults to searched size when no merge/drop occurred
    full_size: Optional[List[int]] = None
    if searched_size is not None:
        full_size = [
            full_m if full_m is not None else searched_size[0],
            full_n if full_n is not None else searched_size[1],
        ]

    # Status classification
    error_obj: Optional[Dict[str, str]] = None
    if killed:
        status = "timeout"
        stop_code = "killed"
        stop_reason = "killed by harness at kill_after_s"
        exact = False
    elif rc is not None and rc != 0:
        status = "error"
        exact = False
        err_msg = _extract_error_msg(clean_err) or f"exit code {rc}"
        error_obj = {"kind": "engine", "msg": err_msg}
    elif lb is None and not k0:
        status = "error"
        exact = False
        error_obj = {"kind": "engine", "msg": "missing result line on stdout"}
    elif k0:
        status = "certified"
        exact = True
        stop_code = stop_code or "k0"
    elif config == "cc_fixed":
        exact = bool(lb is not None and ub is not None and lb > 0 and lb == ub)
        if stop_code == "timeout":
            status = "timeout"
        elif exact:
            status = "certified"
        else:
            status = "bracket"
    elif lb is not None and ub is not None and lb > 0 and lb == ub:
        status = "certified"
        exact = True
    elif stop_code == "timeout":
        status = "timeout"
        exact = False
    else:
        status = "bracket"
        exact = False

    cc_rounds = [cc_rounds_by_w[w] for w in sorted(cc_rounds_by_w)]
    return {
        "status": status,
        "status_engine": status,
        "lb": lb,
        "ub": ub,
        "exact": exact,
        "t_internal": {"t_stop": t_stop, "t_input": None, "t_import": None, "t_call": None},
        "stop_reason": stop_reason,
        "stop_code": stop_code,
        "t_ub": t_ub,
        "t_lb": t_lb,
        "rw_steps": rw_steps,
        "size": {"full": full_size, "searched": searched_size},
        "cc_rounds": cc_rounds,
        "threads_used": {"requested": requested_threads, "rw": rw_threads, "cc": cc_threads},
        "rw_opt": rw_opt,
        "k0": k0,
        "witness": chosen_witness,
        "events": events,
        "parse": "stderr",
        "engine_info": {},
        "error": error_obj,
    }


def _from_json_result(
    obj: Dict[str, Any],
    *,
    rc: Optional[int],
    killed: bool,
    config: str,
    dem_mode: str,
) -> Dict[str, Any]:
    k0 = bool(obj.get("k0", False))
    dmin_raw = obj.get("dmin")
    dmax_raw = obj.get("dmax")
    lb = int(dmin_raw) if dmin_raw is not None else None
    if k0 and dmax_raw == 0:
        ub: Optional[int] = 0
    elif dmax_raw is not None and int(dmax_raw) > 0:
        ub = int(dmax_raw)
    else:
        ub = None

    stop = obj.get("stop") or {}
    stop_code = stop.get("code")
    stop_reason = stop.get("text")

    w_obj = obj.get("witness") or {}
    witness: Optional[Dict[str, Any]] = None
    if w_obj.get("available"):
        base_in = int(w_obj.get("base", 1))
        sup = [int(c) - base_in for c in (w_obj.get("support") or [])]
        wt = int(w_obj.get("weight", len(sup)))
        witness = {
            "weight": wt,
            "support": sup,
            "base": 0,
            "source": w_obj.get("source"),
            "complete": len(sup) == wt and wt > 0,
            "dem_mode": dem_mode,
        }

    events_norm: List[List[Any]] = []
    for ev in obj.get("events") or []:
        if isinstance(ev, list) and len(ev) >= 3:
            ev_ub = int(ev[2]) if (ev[2] is not None and int(ev[2]) > 0) else None
            events_norm.append([ev[0], ev[1], ev_ub])

    thr = obj.get("threads") or {}
    sz = obj.get("size") or {}
    error_obj: Optional[Dict[str, str]] = None
    if killed:
        status = "timeout"
        stop_code = "killed"
        exact = False
    elif rc is not None and rc != 0:
        status = "error"
        exact = False
        error_obj = {"kind": "engine", "msg": f"exit code {rc}"}
    elif k0:
        status = "certified"
        exact = True
    elif config == "cc_fixed":
        exact = bool(lb is not None and ub is not None and lb > 0 and lb == ub)
        status = "timeout" if stop_code == "timeout" else ("certified" if exact else "bracket")
    elif lb is not None and ub is not None and lb > 0 and lb == ub:
        status = "certified"
        exact = True
    elif stop_code == "timeout":
        status = "timeout"
        exact = False
    else:
        status = "bracket"
        exact = False

    return {
        "status": status,
        "status_engine": status,
        "lb": lb,
        "ub": ub,
        "exact": exact,
        "t_internal": {
            "t_stop": obj.get("t_stop"),
            "t_input": obj.get("t_input"),
            "t_import": None,
            "t_call": None,
        },
        "stop_reason": stop_reason,
        "stop_code": stop_code,
        "t_ub": obj.get("t_ub"),
        "t_lb": obj.get("t_lb"),
        "rw_steps": obj.get("rw_steps"),
        "size": {"full": sz.get("full"), "searched": sz.get("searched")},
        "cc_rounds": obj.get("cc_rounds") or [],
        "threads_used": {
            "requested": thr.get("requested"),
            "rw": thr.get("rw"),
            "cc": thr.get("cc"),
        },
        "rw_opt": obj.get("rw_opt"),
        "k0": k0,
        "witness": witness,
        "events": events_norm,
        "parse": "json",
        "engine_info": {
            "search": obj.get("search"),
            "rw_steps_done": obj.get("rw_steps_done"),
            "threads_detail": thr,
        },
        "error": error_obj,
    }


def _extract_error_msg(stderr: str) -> str:
    lines = [ln.strip() for ln in stderr.splitlines() if ln.strip()]
    for i, ln in enumerate(lines):
        if "*** ERROR" in ln and i + 1 < len(lines):
            return lines[i + 1][:240]
        if "unrecognized parameter" in ln:
            return ln[:240]
    for ln in lines:
        if not ln.startswith("#"):
            return ln[:240]
    return ""


def _m4ri_version() -> Optional[str]:
    try:
        p = subprocess.run(["pkg-config", "--modversion", "m4ri"], capture_output=True, text=True, timeout=5)
        if p.returncode == 0 and p.stdout.strip():
            return p.stdout.strip()
    except Exception:
        pass
    return None


class DistM4riBinAdapter(Adapter):
    """Adapter for the native `dist_m4ri` binary (`method=3` default or `method=2` fixed-work scaling)."""

    name = "dist_m4ri"
    interface = "binary"

    def __init__(
        self,
        binary_path: Optional[Path] = None,
        config: str = "default",
        repo_root: Path = REPO_ROOT,
        cache_dir: Optional[Path] = None,
        allow_fallback: bool = True,
    ) -> None:
        if config not in ("default", "cc_fixed"):
            raise ValueError(f"Unsupported dist_m4ri binary config {config!r} (expected 'default' or 'cc_fixed')")
        self.config = config
        self.repo_root = repo_root.resolve()
        self.binary = (binary_path if binary_path is not None else (self.repo_root / "src" / "dist_m4ri")).resolve()
        self.cache_dir = (cache_dir if cache_dir is not None else mi.default_cache_dir(self.repo_root)).resolve()
        self.allow_fallback = allow_fallback
        self._has_json: Optional[bool] = None
        self._desc_cache: Optional[Dict[str, Any]] = None
        self._last_threads: Optional[int] = None

    def probe_json_support(self) -> bool:
        """Probe `<binary> --help json` once; rc == 0 means opt-in json=/stamp= are supported."""
        if self._has_json is not None:
            return self._has_json
        if not self.binary.is_file():
            self._has_json = False
            return False
        try:
            p = subprocess.run(
                [str(self.binary), "--help", "json"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            out = (p.stdout or "") + (p.stderr or "")
            self._has_json = (p.returncode == 0 and "json=" in out)
        except Exception:
            self._has_json = False
        return self._has_json

    def describe(self) -> Dict[str, Any]:
        if self._desc_cache is not None:
            return dict(self._desc_cache)
        ver = None
        if self.binary.is_file():
            try:
                p = subprocess.run([str(self.binary), "--version"], capture_output=True, text=True, timeout=5)
                txt = (p.stdout + p.stderr).strip()
                m = re.search(r"([0-9]+\.[0-9]+\.[0-9]+[^\s]*)", txt)
                ver = m.group(1) if m else txt
            except Exception:
                ver = None
        commit, dirty = git_info(self.repo_root)
        has_json = self.probe_json_support()
        self._desc_cache = {
            "name": self.name,
            "interface": self.interface,
            "config": self.config,
            "version": ver,
            "commit": commit,
            "dirty": dirty,
            "binary": str(self.binary),
            "sha256": mi.sha256_file(self.binary) if self.binary.is_file() else None,
            "deps": {"m4ri": _m4ri_version()},
            "features": {"json": has_json, "stamp": has_json},
        }
        return dict(self._desc_cache)

    def supports(self, inst: Dict[str, Any]) -> bool:
        if not self.binary.is_file():
            return False
        if inst["kind"] not in ("css_sector", "classical", "hl", "dem", "circuit"):
            return False
        if self.config == "cc_fixed":
            d_ref = inst.get("d_ref")
            return d_ref is not None and int(d_ref) >= 2
        return True

    def prepare(self, inst: Dict[str, Any]) -> Dict[str, Any]:
        files = mi.resolve_instance_files(
            inst, repo_root=self.repo_root, cache_dir=self.cache_dir, verify_sha=False
        )
        if inst["kind"] != "circuit":
            return {"prep_s": 0.0, "files": {k: str(v) for k, v in files.items()}, "dem_mode": "canonical"}

        # Native binary cannot read .stim directly; convert to canonical DEM in cache and record prep_s
        recipe = inst.get("recipe")
        if not recipe:
            raise ValueError(f"Circuit instance {inst['id']} lacks a recipe to generate canonical DEM")
        src_sha = recipe["source_stim"]["sha256"]
        prep_dir = self.cache_dir / "prep_dem" / src_sha[:2]
        dem_path = prep_dir / f"{src_sha}_canonical.dem"
        meta_path = prep_dir / f"{src_sha}_canonical.meta.json"
        if dem_path.is_file() and meta_path.is_file():
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            if mi.sha256_file(dem_path) == meta.get("sha256"):
                return {
                    "prep_s": float(meta.get("gen_s", 0.0)),
                    "files": {"dem": str(dem_path)},
                    "dem_mode": "canonical",
                    "dem_path": str(dem_path),
                }
        prep_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.monotonic()
        data, meta = mi.generate_canonical_dem(recipe, repo_root=self.repo_root)
        gen_s = round(time.monotonic() - t0, 4)
        meta["gen_s"] = gen_s
        tmp = dem_path.with_suffix(".dem.tmp")
        tmp.write_bytes(data)
        os.replace(tmp, dem_path)
        meta_tmp = meta_path.with_suffix(".json.tmp")
        meta_tmp.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        os.replace(meta_tmp, meta_path)
        return {
            "prep_s": gen_s,
            "files": {"dem": str(dem_path)},
            "dem_mode": "canonical",
            "dem_path": str(dem_path),
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
        files = prep["files"]
        kind = inst["kind"]
        if kind == "css_sector":
            in_args = [f"finH={files['H']}", f"finG={files['G']}"]
        elif kind == "classical":
            in_args = [f"finH={files['H']}", "classical=1"]
        elif kind == "hl":
            in_args = [f"finH={files['H']}", f"finL={files['L']}"]
        elif kind in ("dem", "circuit"):
            in_args = [f"fdem={files['dem']}"]
        else:
            raise ValueError(f"Unsupported kind {kind!r}")

        t_str = f"{budget_s:g}"
        if self.config == "cc_fixed":
            wmax = int(inst["d_ref"]) - 1
            extra = ["method=2", f"wmax={wmax}", f"threads={threads}", f"timeout={t_str}", f"seed={seed}"]
        else:
            extra = ["method=3", f"threads={threads}", f"timeout={t_str}", f"seed={seed}"]

        # DBG_SUMMARY(1) | DBG_PROGRESS(2) | DBG_CODEWORDS(32); add DBG_LARGE(2048) if d may exceed 25.
        # Kept active even with json=1 stamp=1 so hard-killed runs still emit verifiable intermediate witnesses.
        d_ref = inst.get("d_ref")
        br = inst.get("bracket_ref")
        may_exceed_25 = (
            (d_ref is not None and int(d_ref) > 25)
            or (br is not None and int(br[1]) > 25)
            or (d_ref is None and br is None)
        )
        dbg = (35 | 2048) if may_exceed_25 else 35
        if self.probe_json_support():
            extra.extend(["json=1", "stamp=1", f"debug={dbg}"])
        else:
            extra.append(f"debug={dbg}")

        self._last_threads = threads
        argv = [str(self.binary)] + in_args + extra
        return {"argv": argv, "argv_extra": extra, "env": {}, "cwd": str(self.repo_root)}

    def parse(
        self,
        rc: Optional[int],
        stdout: str,
        stderr: str,
        *,
        killed: bool = False,
        prep: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        if rc == 255 and "unrecognized parameter" in (stderr or "") and self.probe_json_support():
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
                "parse": "error",
                "engine_info": {},
                "error": {
                    "kind": "harness",
                    "msg": "binary rejected json=/stamp= despite --help json returning 0",
                },
            }
        dem_mode = (prep or {}).get("dem_mode", "canonical")
        return parse_binary_stderr(
            stdout,
            stderr,
            rc=rc,
            killed=killed,
            config=self.config,
            dem_mode=dem_mode,
            requested_threads=self._last_threads,
        )
