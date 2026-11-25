#!/usr/bin/env python3
"""Independent GF(2) witness verification and lower-bound cross-audit for dist-m4ri benchmarks.

Checks every run record (`dm4ri-bench-run/1`) outside the timed region:
  1. Independent numpy GF(2) check of the reported witness support:
     - classical:  H e = 0 (mod 2), e != 0
     - css_sector: H e = 0 (mod 2), e not in rowspan(G)
     - hl:         H e = 0 (mod 2), L e != 0 (mod 2)
     - dem/circ:   H_dem e = 0 (mod 2), L_dem e != 0 (mod 2)
  2. Reference and self-consistency checks:
     - lb <= ub (when ub > 0), lb <= verified witness weight, lb <= d_ref (or bracket_ref[1])
     - certified distance equals d_ref (or lies inside bracket_ref)
     - verified witness weight >= d_ref (or bracket_ref[0]); a lighter valid witness flags a conflict
  3. Cross-run lower-bound audit (`audit_records`):
     - across all runs of an instance, every lb <= min verified witness weight

A failed check sets `result.status = "wrong"`, records the details in `verify.violations`, and marks
`verify.lb_audit = "violation"` where a lower bound is implicated.

CLI usage:
  python3 tools/bench/verify.py --manifest benchmark/manifest/certify_v1.json --session <dir>
  python3 tools/bench/verify.py --manifest benchmark/manifest/certify_v1.json --session <dir> --audit
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST_DIR = REPO_ROOT / "benchmark" / "manifest"
if str(MANIFEST_DIR) not in sys.path:
    sys.path.insert(0, str(MANIFEST_DIR))
import make_inputs as mi  # noqa: E402


# ---------------------------------------------------------------------------
# MatrixMarket & DEM readers (independent of dist_m4ri C binary)
# ---------------------------------------------------------------------------

def read_mtx_gf2(path: Path) -> np.ndarray:
    """Read a 1-based MatrixMarket coordinate file into a uint8 GF(2) matrix (duplicates XORed)."""
    nrows: Optional[int] = None
    ncols: Optional[int] = None
    mat: Optional[np.ndarray] = None
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"{path}: missing %%MatrixMarket header")
        hdr_lo = header.lower()
        is_pattern = "pattern" in hdr_lo
        is_sym = "symmetric" in hdr_lo
        for line in f:
            s = line.strip()
            if not s or s.startswith("%"):
                continue
            parts = s.split()
            if nrows is None:
                nrows, ncols = int(parts[0]), int(parts[1])
                mat = np.zeros((nrows, ncols), dtype=np.uint8)
                continue
            r = int(parts[0]) - 1
            c = int(parts[1]) - 1
            val = 1 if is_pattern or len(parts) < 3 else (int(float(parts[2])) & 1)
            if val and mat is not None:
                mat[r, c] ^= 1
                if is_sym and r != c:
                    mat[c, r] ^= 1
    if mat is None:
        raise ValueError(f"{path}: empty MatrixMarket file")
    return mat


def _pack_rows_u64(mat: np.ndarray) -> np.ndarray:
    """Pack a 2D uint8 GF(2) matrix into uint64 bit-words per row."""
    arr = np.asarray(mat, dtype=np.uint8) & 1
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    m, n = arr.shape
    if n == 0:
        return np.zeros((m, 0), dtype=np.uint64)
    packed_bytes = np.packbits(arr, axis=1, bitorder="little")
    rem = packed_bytes.shape[1] % 8
    if rem:
        pad = np.zeros((m, 8 - rem), dtype=np.uint8)
        packed_bytes = np.hstack([packed_bytes, pad])
    return np.ascontiguousarray(packed_bytes).view(np.uint64)


def _leading_bit(row_u64: np.ndarray) -> int:
    """Return index of the lowest set bit in a uint64 word array, or -1 if zero."""
    nz = np.flatnonzero(row_u64)
    if nz.size == 0:
        return -1
    w_idx = int(nz[0])
    word = int(row_u64[w_idx])
    lsb = word & -word
    return (w_idx << 6) + (lsb.bit_length() - 1)


class GF2RowBasis:
    """Precomputed GF(2) row-echelon basis for fast membership testing (e in rowspan(G))."""

    def __init__(self, mat: np.ndarray) -> None:
        self.ncols = int(mat.shape[1])
        self.pivots: Dict[int, np.ndarray] = {}
        packed = _pack_rows_u64(mat)
        for i in range(packed.shape[0]):
            self._insert(packed[i].copy())

    def _insert(self, row: np.ndarray) -> bool:
        while True:
            bit = _leading_bit(row)
            if bit < 0:
                return False
            piv = self.pivots.get(bit)
            if piv is None:
                self.pivots[bit] = row
                return True
            np.bitwise_xor(row, piv, out=row)

    @property
    def rank(self) -> int:
        return len(self.pivots)

    def contains_support(self, support: Sequence[int]) -> bool:
        """Return True iff the binary vector with 1s at `support` lies in the row space."""
        vec = np.zeros((1, self.ncols), dtype=np.uint8)
        if support:
            vec[0, list(support)] = 1
        row = _pack_rows_u64(vec)[0].copy()
        while True:
            bit = _leading_bit(row)
            if bit < 0:
                return True
            piv = self.pivots.get(bit)
            if piv is None:
                return False
            np.bitwise_xor(row, piv, out=row)


def _dem_inst_parts(line: str) -> Tuple[str, str, str]:
    """Parse one DEM line into (instruction_name, arg_str, target_str), stripping tags and comments."""
    s = line.strip()
    if not s or s[0] == "#":
        return "", "", ""
    if s[0] == "}":
        return "}", "", ""
    i = 0
    while i < len(s) and (s[i].isalnum() or s[i] == "_"):
        i += 1
    name, j = s[:i], i
    if j < len(s) and s[j] == "[":
        j += 1
        while j < len(s) and s[j] != "]":
            j += 2 if s[j] == "\\" else 1
        j = min(j + 1, len(s))
    rest = s[j:].split("#", 1)[0].strip()
    if rest.startswith("("):
        k = rest.find(")")
        if k < 0:
            return name, rest[1:].strip(), ""
        return name, rest[1:k].strip(), rest[k + 1:].strip()
    return name, "", rest


def parse_dem_columns(dem_text: str) -> Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], int, int]:
    """Parse a DEM (with repeat/shift_detectors/tags/^) into (cols, num_detectors, num_observables).

    Matches the dist_m4ri C binary's column numbering: every error(...) line in the unrolled DEM
    corresponds to one column (0-based), with `^`-separated components summed mod 2.
    """
    lines = dem_text.splitlines()
    cols: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []
    max_d = -1
    max_o = -1

    def run_block(idx: int, shift: int) -> Tuple[int, int]:
        nonlocal max_d, max_o
        while idx < len(lines):
            name, _args, rest = _dem_inst_parts(lines[idx])
            idx += 1
            if not name:
                continue
            if name == "}":
                return idx, shift
            if name == "repeat":
                count = int(rest.split()[0])
                body = idx
                end = idx
                for _ in range(count):
                    end, shift = run_block(body, shift)
                if count <= 0:
                    depth = 1
                    while end < len(lines) and depth > 0:
                        nm = _dem_inst_parts(lines[end])[0]
                        end += 1
                        depth += (nm == "repeat") - (nm == "}")
                idx = end
            elif name == "shift_detectors":
                shift += int(rest.split()[-1])
            elif name == "detector":
                for tok in rest.split():
                    if tok.startswith("D"):
                        max_d = max(max_d, int(tok[1:]) + shift)
            elif name == "logical_observable":
                for tok in rest.split():
                    if tok.startswith("L"):
                        max_o = max(max_o, int(tok[1:]))
            elif name == "error":
                dets: set[int] = set()
                obs: set[int] = set()
                for tok in rest.split():
                    if tok == "^":
                        continue
                    if tok[0] == "D":
                        d = int(tok[1:]) + shift
                        max_d = max(max_d, d)
                        dets ^= {d}
                    elif tok[0] == "L":
                        o = int(tok[1:])
                        max_o = max(max_o, o)
                        obs ^= {o}
                cols.append((tuple(sorted(dets)), tuple(sorted(obs))))
        return idx, shift

    run_block(0, 0)
    return cols, max_d + 1, max_o + 1


# ---------------------------------------------------------------------------
# Witness verifiers
# ---------------------------------------------------------------------------

def verify_matrix_witness_numpy(
    H: np.ndarray,
    support: Sequence[int],
    *,
    L: Optional[np.ndarray] = None,
    G: Optional[np.ndarray] = None,
    G_basis: Optional[GF2RowBasis] = None,
    classical: bool = False,
) -> Dict[str, Any]:
    """Verify a 0-based codeword support over GF(2) using numpy."""
    H_u8 = np.asarray(H, dtype=np.uint8) & 1
    ncols = int(H_u8.shape[1])
    cols = [int(c) for c in support]
    if not cols:
        return {"ok": False, "reason": "empty witness", "weight": 0}
    if len(set(cols)) != len(cols):
        return {"ok": False, "reason": "duplicate column indices", "weight": len(cols)}
    if any(c < 0 or c >= ncols for c in cols):
        return {"ok": False, "reason": f"column index out of range [0, {ncols})", "weight": len(cols)}

    h_syn = np.bitwise_xor.reduce(H_u8[:, cols], axis=1) if H_u8.shape[0] > 0 else np.zeros(0, dtype=np.uint8)
    h_wt = int(np.count_nonzero(h_syn))
    if h_wt != 0:
        return {
            "ok": False,
            "reason": f"H e != 0 (syndrome weight {h_wt})",
            "weight": len(cols),
            "H_syndrome_weight": h_wt,
        }

    detail: Dict[str, Any] = {"ok": True, "weight": len(cols), "H_syndrome_weight": 0}
    if classical:
        return detail

    checked_nontrivial = False
    if L is not None:
        L_u8 = np.asarray(L, dtype=np.uint8) & 1
        l_syn = np.bitwise_xor.reduce(L_u8[:, cols], axis=1) if L_u8.shape[0] > 0 else np.zeros(0, dtype=np.uint8)
        l_wt = int(np.count_nonzero(l_syn))
        detail["L_syndrome_weight"] = l_wt
        checked_nontrivial = True
        if l_wt == 0:
            detail["ok"] = False
            detail["reason"] = "L e == 0 (mod 2)"
            return detail

    if G_basis is not None or G is not None:
        basis = G_basis if G_basis is not None else GF2RowBasis(G)  # type: ignore[arg-type]
        in_rowspan = basis.contains_support(cols)
        detail["independent_of_G"] = not in_rowspan
        checked_nontrivial = True
        if in_rowspan:
            detail["ok"] = False
            detail["reason"] = "e lies in rowspan(G)"
            return detail

    if not checked_nontrivial:
        return {"ok": False, "reason": "neither L nor G provided for non-triviality check", "weight": len(cols)}
    return detail


def verify_dem_witness_numpy(
    cols: Sequence[Tuple[Tuple[int, ...], Tuple[int, ...]]],
    num_detectors: int,
    num_observables: int,
    support: Sequence[int],
) -> Dict[str, Any]:
    """Verify a 0-based DEM error-column support over GF(2) using numpy syndrome vectors."""
    ncols = len(cols)
    idx_list = [int(c) for c in support]
    if not idx_list:
        return {"ok": False, "reason": "empty witness", "weight": 0}
    if len(set(idx_list)) != len(idx_list):
        return {"ok": False, "reason": "duplicate column indices", "weight": len(idx_list)}
    if any(c < 0 or c >= ncols for c in idx_list):
        return {"ok": False, "reason": f"column index out of range [0, {ncols})", "weight": len(idx_list)}

    h_syn = np.zeros(max(0, num_detectors), dtype=np.uint8)
    l_syn = np.zeros(max(0, num_observables), dtype=np.uint8)
    for c in idx_list:
        dets, obs = cols[c]
        if dets:
            h_syn[list(dets)] ^= 1
        if obs:
            l_syn[list(obs)] ^= 1

    h_wt = int(np.count_nonzero(h_syn))
    l_wt = int(np.count_nonzero(l_syn))
    if h_wt != 0:
        return {
            "ok": False,
            "reason": f"H_dem e != 0 (detector syndrome weight {h_wt})",
            "weight": len(idx_list),
            "H_syndrome_weight": h_wt,
            "L_syndrome_weight": l_wt,
        }
    if l_wt == 0:
        return {
            "ok": False,
            "reason": "L_dem e == 0 (no logical observable flipped)",
            "weight": len(idx_list),
            "H_syndrome_weight": 0,
            "L_syndrome_weight": 0,
        }
    return {
        "ok": True,
        "weight": len(idx_list),
        "H_syndrome_weight": 0,
        "L_syndrome_weight": l_wt,
    }


# ---------------------------------------------------------------------------
# Verifier cache & record-level check
# ---------------------------------------------------------------------------

class VerifierContext:
    """Caches loaded matrices, GF(2) row bases, and parsed DEM columns across runs."""

    def __init__(self, repo_root: Path = REPO_ROOT, cache_dir: Optional[Path] = None) -> None:
        self.repo_root = repo_root.resolve()
        self.cache_dir = (cache_dir if cache_dir is not None else mi.default_cache_dir(self.repo_root)).resolve()
        self._mtx_cache: Dict[Path, np.ndarray] = {}
        self._basis_cache: Dict[Path, GF2RowBasis] = {}
        self._dem_cache: Dict[str, Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], int, int]] = {}

    def _get_mtx(self, path: Path) -> np.ndarray:
        mat = self._mtx_cache.get(path)
        if mat is None:
            mat = read_mtx_gf2(path)
            self._mtx_cache[path] = mat
        return mat

    def _get_basis(self, path: Path) -> GF2RowBasis:
        basis = self._basis_cache.get(path)
        if basis is None:
            basis = GF2RowBasis(self._get_mtx(path))
            self._basis_cache[path] = basis
        return basis

    def _get_dem_file(self, path: Path) -> Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], int, int]:
        key = f"file:{path.resolve()}"
        val = self._dem_cache.get(key)
        if val is None:
            val = parse_dem_columns(path.read_text(encoding="utf-8"))
            # Avoid retaining huge DEM column lists in the harness process (preserves wait4 max_rss floor).
            if len(val[0]) <= 50000:
                self._dem_cache[key] = val
        return val

    def _get_circuit_dem(
        self,
        inst: Dict[str, Any],
        dem_mode: str,
        dem_path_hint: Optional[str] = None,
    ) -> Tuple[List[Tuple[Tuple[int, ...], Tuple[int, ...]]], int, int]:
        if dem_path_hint and Path(dem_path_hint).is_file():
            return self._get_dem_file(Path(dem_path_hint))
        key = f"circ:{inst['id']}:{dem_mode}"
        val = self._dem_cache.get(key)
        if val is not None:
            return val
        recipe = inst.get("recipe")
        if not recipe:
            raise ValueError(f"Circuit instance {inst['id']} has no recipe for DEM verification")
        if dem_mode == "stim_auto":
            # Matches dist_m4ri.py _dem_from_input -> _circuit_dem (decompose_errors=True fallback False)
            dm = mi._import_dist_m4ri(self.repo_root)
            stim = dm._get_stim()
            src_path = self.repo_root / recipe["source_stim"]["path"]
            circ = stim.Circuit.from_file(str(src_path))
            circ, _ = dm.remove_empty_detectors(circ)
            if recipe.get("rounds") is not None:
                circ = dm.set_circuit_rounds(circ, int(recipe["rounds"]))
            if recipe.get("add_noise") and not dm.has_noise(circ):
                circ = dm.add_noise(circ, p=float(recipe["add_noise"]["p"]))
            dem, _, _ = dm._circuit_dem(circ)
            text = str(dem.flattened())
            val = parse_dem_columns(text)
        else:
            data, _ = mi.generate_canonical_dem(recipe, repo_root=self.repo_root)
            val = parse_dem_columns(data.decode("utf-8"))
        if len(val[0]) <= 50000:
            self._dem_cache[key] = val
        return val

    def check_witness(
        self,
        inst: Dict[str, Any],
        witness: Dict[str, Any],
        dem_path_hint: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run the independent numpy GF(2) check for a witness dict on `inst`."""
        support_raw = witness.get("support") or []
        base = int(witness.get("base", 0))
        support = [int(c) - base for c in support_raw]
        kind = inst["kind"]
        files = mi.resolve_instance_files(inst, repo_root=self.repo_root, cache_dir=self.cache_dir, verify_sha=False)

        if kind == "classical":
            return verify_matrix_witness_numpy(self._get_mtx(files["H"]), support, classical=True)
        if kind == "css_sector":
            return verify_matrix_witness_numpy(
                self._get_mtx(files["H"]), support, G_basis=self._get_basis(files["G"])
            )
        if kind == "hl":
            return verify_matrix_witness_numpy(self._get_mtx(files["H"]), support, L=self._get_mtx(files["L"]))
        if kind == "dem":
            cols, m, k = self._get_dem_file(files["dem"])
            return verify_dem_witness_numpy(cols, m, k, support)
        if kind == "circuit":
            dem_mode = str(witness.get("dem_mode") or "canonical")
            cols, m, k = self._get_circuit_dem(inst, dem_mode, dem_path_hint=dem_path_hint)
            return verify_dem_witness_numpy(cols, m, k, support)
        raise ValueError(f"Unsupported instance kind {kind!r}")

    def verify_record(
        self,
        record: Dict[str, Any],
        inst: Dict[str, Any],
        dem_path_hint: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Verify a run record in place, populating `record['verify']` and setting `status='wrong'` on violation."""
        res = record["result"]
        status_engine = res.get("status_engine") or res.get("status")
        res["status_engine"] = status_engine
        lb = res.get("lb")
        ub = res.get("ub")
        exact = bool(res.get("exact"))
        witness = res.get("witness")

        d_ref = inst.get("d_ref")
        bracket_ref = inst.get("bracket_ref")

        violations: List[str] = []
        witness_weight: Optional[int] = None
        witness_verified: Optional[bool] = None
        witness_detail: Optional[Dict[str, Any]] = None
        ub_verified: Optional[bool] = None
        unverified_reason: Optional[str] = None
        lb_audit = "ok" if lb is not None else "n/a"
        matches_ref: Optional[bool] = None

        # 1. Witness GF(2) verification
        if witness is not None:
            w_declared = witness.get("weight")
            w_complete = bool(witness.get("complete", True))
            w_support = witness.get("support") or []
            if not w_complete or not w_support:
                unverified_reason = (
                    "witness support truncated or incomplete"
                    if not w_complete
                    else "empty witness support"
                )
                ub_verified = False if (ub is not None and ub > 0) else None
            else:
                detail = self.check_witness(inst, witness, dem_path_hint=dem_path_hint)
                witness_detail = detail
                witness_weight = int(detail.get("weight", len(w_support)))
                witness_verified = bool(detail.get("ok"))
                if not witness_verified:
                    violations.append(f"witness failed GF(2) check: {detail.get('reason', 'unknown')}")
                elif w_declared is not None and int(w_declared) != witness_weight:
                    violations.append(
                        f"witness declared weight {w_declared} != support length {witness_weight}"
                    )
                    witness_verified = False
                if ub is not None and ub > 0:
                    ub_verified = bool(witness_verified and witness_weight == ub)
                    if witness_verified and witness_weight != ub:
                        unverified_reason = f"verified witness has weight {witness_weight} != reported ub {ub}"
        elif ub is not None and ub > 0:
            ub_verified = False
            unverified_reason = str(
                (res.get("engine_info") or {}).get("witness_unavailable_reason")
                or "no witness reported by engine"
            )

        # 2. Bound self-consistency & witness vs lb
        if lb is not None and ub is not None and ub > 0 and lb > ub:
            violations.append(f"reported lb ({lb}) > reported ub ({ub})")
            lb_audit = "violation"
        if witness_verified is True and witness_weight is not None and lb is not None and lb > witness_weight:
            violations.append(f"reported lb ({lb}) > verified witness weight ({witness_weight})")
            lb_audit = "violation"

        # 3. Reference check (d_ref / bracket_ref)
        if status_engine != "error":
            if d_ref is not None:
                if lb is not None and lb > d_ref:
                    violations.append(f"reported lb ({lb}) > d_ref ({d_ref})")
                    lb_audit = "violation"
                if witness_verified is True and witness_weight is not None and witness_weight < d_ref:
                    violations.append(
                        f"verified witness weight ({witness_weight}) < d_ref ({d_ref}); check reference"
                    )
                if exact or status_engine == "certified":
                    matches_ref = (lb == d_ref and ub == d_ref)
                    if not matches_ref:
                        violations.append(f"certified [{lb}, {ub}] != d_ref ({d_ref})")
                        lb_audit = "violation" if (lb is not None and lb > d_ref) else lb_audit
            elif bracket_ref is not None:
                ref_lo, ref_hi = int(bracket_ref[0]), int(bracket_ref[1])
                if lb is not None and lb > ref_hi:
                    violations.append(f"reported lb ({lb}) > bracket_ref[1] ({ref_hi})")
                    lb_audit = "violation"
                if witness_verified is True and witness_weight is not None and witness_weight < ref_lo:
                    violations.append(
                        f"verified witness weight ({witness_weight}) < bracket_ref[0] ({ref_lo}); check reference"
                    )
                if exact or status_engine == "certified":
                    matches_ref = (lb is not None and ref_lo <= lb <= ref_hi and lb == ub)
                    if not matches_ref:
                        violations.append(
                            f"certified [{lb}, {ub}] outside bracket_ref [{ref_lo}, {ref_hi}]"
                        )
                        if lb is not None and lb > ref_hi:
                            lb_audit = "violation"

        verify_block = {
            "witness_weight": witness_weight,
            "witness_verified_numpy": witness_verified,
            "witness_check_detail": witness_detail,
            "ub_verified": ub_verified,
            "lb_audit": lb_audit,
            "matches_ref": matches_ref,
            "violations": violations,
            "unverified_reason": unverified_reason,
        }
        record["verify"] = verify_block
        if violations:
            res["status"] = "wrong"
        else:
            res["status"] = status_engine
        return record


def deduplicate_records(records: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Keep the last record for each run_id, preserving last-seen order."""
    by_id: Dict[str, Dict[str, Any]] = {}
    for rec in records:
        rid = rec.get("run_id")
        if rid is not None:
            by_id[str(rid)] = rec
    return list(by_id.values())


def audit_records(
    records: Sequence[Dict[str, Any]],
    manifest: Dict[str, Any],
) -> Dict[str, Any]:
    """Cross-run integrity audit across all deduplicated records (C4 / Table C).

    Checks:
      - wrong runs (status == 'wrong')
      - cross-run LB conflicts: any run with lb > min verified witness weight of the same instance,
        or max lb > d_ref / bracket_ref[1]
      - certified runs without a verified witness
      - contaminated runs
    """
    dedup = deduplicate_records(records)
    by_inst_id = {inst["id"]: inst for inst in manifest.get("instances", [])}
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for rec in dedup:
        iid = rec["instance"]["id"]
        grouped.setdefault(iid, []).append(rec)

    wrong_runs: List[Dict[str, Any]] = []
    lb_conflicts: List[Dict[str, Any]] = []
    unverified_certified: List[Dict[str, Any]] = []
    contaminated_runs: List[Dict[str, Any]] = []

    for iid, recs in sorted(grouped.items()):
        inst = by_inst_id.get(iid, {})
        d_ref = inst.get("d_ref")
        bracket_ref = inst.get("bracket_ref")
        verified_weights = [
            int(r["verify"]["witness_weight"])
            for r in recs
            if (r.get("verify") or {}).get("witness_verified_numpy") is True
            and (r.get("verify") or {}).get("witness_weight") is not None
        ]
        min_ver_w = min(verified_weights) if verified_weights else None

        for r in recs:
            res = r.get("result") or {}
            ver = r.get("verify") or {}
            env = r.get("env") or {}
            status = res.get("status")
            lb = res.get("lb")

            if status == "wrong":
                wrong_runs.append({
                    "run_id": r["run_id"],
                    "instance": iid,
                    "engine": f"{r['engine']['name']}:{r['engine']['config']}@p{r['options']['threads']}",
                    "seed": r["options"]["seed"],
                    "lb": lb,
                    "ub": res.get("ub"),
                    "violations": ver.get("violations") or [],
                })
            if res.get("status_engine") != "error" and lb is not None:
                conflict_reasons: List[str] = []
                if min_ver_w is not None and lb > min_ver_w:
                    conflict_reasons.append(f"lb ({lb}) > cross-run min verified witness ({min_ver_w})")
                if d_ref is not None and lb > d_ref:
                    conflict_reasons.append(f"lb ({lb}) > d_ref ({d_ref})")
                elif bracket_ref is not None and lb > int(bracket_ref[1]):
                    conflict_reasons.append(f"lb ({lb}) > bracket_ref[1] ({bracket_ref[1]})")
                if conflict_reasons:
                    lb_conflicts.append({
                        "run_id": r["run_id"],
                        "instance": iid,
                        "engine": f"{r['engine']['name']}:{r['engine']['config']}@p{r['options']['threads']}",
                        "seed": r["options"]["seed"],
                        "lb": lb,
                        "min_verified_witness": min_ver_w,
                        "reasons": conflict_reasons,
                    })
            ub_ok = ver.get("ub_verified") if "ub_verified" in ver else ver.get("witness_verified_numpy")
            if status == "certified" and not res.get("k0") and ub_ok is not True:
                unverified_certified.append({
                    "run_id": r["run_id"],
                    "instance": iid,
                    "engine": f"{r['engine']['name']}:{r['engine']['config']}@p{r['options']['threads']}",
                    "seed": r["options"]["seed"],
                    "ub": res.get("ub"),
                    "reason": ver.get("unverified_reason") or "ub witness not verified",
                })
            if env.get("contaminated"):
                contaminated_runs.append({
                    "run_id": r["run_id"],
                    "instance": iid,
                    "engine": f"{r['engine']['name']}:{r['engine']['config']}@p{r['options']['threads']}",
                    "seed": r["options"]["seed"],
                    "ext_cpu_s": env.get("ext_cpu_s"),
                    "wall_sec": res.get("wall_sec"),
                })

    return {
        "total_runs": len(dedup),
        "wrong_runs": wrong_runs,
        "lb_conflicts": lb_conflicts,
        "unverified_certified": unverified_certified,
        "contaminated_runs": contaminated_runs,
        "ok": not wrong_runs and not lb_conflicts,
    }


def _load_session_records(session_dir: Path) -> List[Dict[str, Any]]:
    res_file = session_dir / "results.jsonl"
    out: List[Dict[str, Any]] = []
    if res_file.is_file():
        for lineno, line in enumerate(res_file.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            try:
                out.append(json.loads(line))
            except ValueError as exc:
                print(
                    f"Warning: skipping corrupt JSON line {lineno} in {res_file}: {exc}",
                    file=sys.stderr,
                )
        return out
    runs_dir = session_dir / "runs"
    if runs_dir.is_dir():
        for rpath in sorted(runs_dir.glob("*/run.json")):
            try:
                out.append(json.loads(rpath.read_text(encoding="utf-8")))
            except (OSError, ValueError) as exc:
                print(f"Warning: skipping unreadable {rpath}: {exc}", file=sys.stderr)
    return out


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Verify witnesses and audit lower bounds for benchmark sessions.")
    p.add_argument("--manifest", type=Path, default=MANIFEST_DIR / "certify_v1.json")
    p.add_argument("--session", type=Path, action="append", required=True, help="Session directory (repeatable).")
    p.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    p.add_argument("--cache-dir", type=Path, default=None)
    p.add_argument("--reverify", action="store_true", help="Re-run GF(2) verification on all records.")
    p.add_argument("--audit", action="store_true", help="Run cross-run lower-bound and integrity audit.")
    p.add_argument("--json-out", type=Path, default=None, help="Write audit report JSON.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    manifest, _ = mi.load_manifest(args.manifest.resolve())
    by_id = {inst["id"]: inst for inst in manifest["instances"]}
    ctx = VerifierContext(repo_root=args.repo_root, cache_dir=args.cache_dir)

    all_records: List[Dict[str, Any]] = []
    for sdir in args.session:
        recs = _load_session_records(sdir.resolve())
        if args.reverify:
            for rec in recs:
                inst = by_id.get(rec["instance"]["id"])
                if inst is not None:
                    ctx.verify_record(rec, inst)
        all_records.extend(recs)

    report = audit_records(all_records, manifest)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2) + "\n")

    print(
        f"Audit ({report['total_runs']} runs): "
        f"wrong={len(report['wrong_runs'])}, "
        f"lb_conflicts={len(report['lb_conflicts'])}, "
        f"unverified_certified={len(report['unverified_certified'])}, "
        f"contaminated={len(report['contaminated_runs'])}"
    )
    for w in report["wrong_runs"]:
        print(f"  WRONG {w['instance']} ({w['engine']} s={w['seed']}): {w['violations']}")
    for c in report["lb_conflicts"]:
        print(f"  LB_CONFLICT {c['instance']} ({c['engine']} s={c['seed']}): {c['reasons']}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
