#!/usr/bin/env python3
"""Deterministic input generation and sha256 verification for dist-m4ri benchmark manifests.

Generates canonical Detector Error Models (DEMs) from tracked Stim circuits into a content-addressed cache
outside git (default: <repo>/scratch/bench_cache/), verifies every tracked and generated file against the
sha256 recorded in the manifest, and refuses to proceed on any mismatch.

Canonical DEM format (matching scratch/bench_0.14.0/run_bench.py:228-263):
  - remove_empty_detectors(circuit)
  - optional set_circuit_rounds(circuit, rounds) when recipe['rounds'] is not None
  - optional add_noise(circuit, p=recipe['add_noise']['p']) when circuit has no noise
  - circuit.detector_error_model(decompose_errors=False, approximate_disjoint_errors=...)
  - flattened error columns via dist_m4ri._dem_error_columns(dem, with_prob=True)
  - one `error(<p>) D... L...` line per non-empty mechanism (sorted detector and observable indices,
    decomposed `^` components summed mod 2, no detector coordinate lines), plus trailing sentinel
    `detector D<m-1>` / `logical_observable L<k-1>` only when the highest index does not appear in any error.

Usage:
  python3 benchmark/manifest/make_inputs.py --manifest benchmark/manifest/certify_v1.json
  python3 benchmark/manifest/make_inputs.py --only-small --check-overlap
  python3 benchmark/manifest/make_inputs.py --verify-only --tier T0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_MANIFEST = Path(__file__).resolve().parent / "certify_v1.json"


def default_cache_dir(repo_root: Path = REPO_ROOT) -> Path:
    """Default content-addressed cache directory outside git."""
    env = os.environ.get("DM4RI_BENCH_CACHE")
    if env:
        return Path(env)
    return repo_root / "scratch" / "bench_cache"


def sha256_bytes(data: bytes) -> str:
    """Return hex sha256 digest of bytes."""
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    """Return hex sha256 digest of a file, streaming in 1 MiB chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest(manifest_path: Path = DEFAULT_MANIFEST) -> Tuple[Dict[str, Any], str]:
    """Load a manifest JSON file and return (manifest_dict, manifest_sha256)."""
    raw = manifest_path.read_bytes()
    data = json.loads(raw.decode("utf-8"))
    if data.get("schema") != "dm4ri-bench-manifest/1":
        raise ValueError(
            f"Unsupported manifest schema {data.get('schema')!r} in {manifest_path} "
            f"(expected 'dm4ri-bench-manifest/1')"
        )
    return data, sha256_bytes(raw)


def _import_dist_m4ri(repo_root: Path = REPO_ROOT) -> Any:
    """Import dist_m4ri from repo_root without mutating caller sys.path permanently."""
    root_str = str(repo_root.resolve())
    if root_str not in sys.path:
        sys.path.insert(0, root_str)
    import dist_m4ri as dm  # noqa: E402

    return dm


def canonical_dem_bytes(dem: Any, dm: Any) -> Tuple[bytes, Dict[str, int]]:
    """Serialize a stim.DetectorErrorModel into the canonical DEM form (run_bench.py:228-263)."""
    cols = dm._dem_error_columns(dem, with_prob=True)
    m = int(dem.num_detectors)
    k = int(dem.num_observables)
    lines: List[str] = []
    max_d, max_o = -1, -1
    for det_set, obs_set, prob in cols:
        if not det_set and not obs_set:
            continue
        dets = sorted(det_set)
        obs = sorted(obs_set)
        if dets:
            max_d = max(max_d, dets[-1])
        if obs:
            max_o = max(max_o, obs[-1])
        lines.append(
            f"error({prob})"
            + "".join(f" D{d}" for d in dets)
            + "".join(f" L{o}" for o in obs)
            + "\n"
        )
    if max_d < m - 1:
        lines.append(f"detector D{m - 1}\n")
    if max_o < k - 1:
        lines.append(f"logical_observable L{k - 1}\n")
    data = "".join(lines).encode("utf-8")
    meta = {"detectors": m, "observables": k, "errors": len(lines) - int(max_d < m - 1) - int(max_o < k - 1)}
    return data, meta


def generate_canonical_dem(
    recipe: Dict[str, Any],
    repo_root: Path = REPO_ROOT,
    allow_stim_mismatch: bool = False,
) -> Tuple[bytes, Dict[str, Any]]:
    """Generate canonical DEM bytes deterministically from a manifest recipe.

    Verifies the source .stim sha256 and the pinned stim version before running Stim.
    """
    if recipe.get("type") != "canonical_dem_v1":
        raise ValueError(f"Unsupported recipe type: {recipe.get('type')!r}")

    src_spec = recipe["source_stim"]
    src_path = repo_root / src_spec["path"]
    if not src_path.is_file():
        raise FileNotFoundError(f"Source Stim circuit not found: {src_path}")
    actual_src_sha = sha256_file(src_path)
    expected_src_sha = src_spec.get("sha256")
    if expected_src_sha and actual_src_sha != expected_src_sha:
        raise ValueError(
            f"Source .stim sha256 mismatch for {src_spec['path']}: "
            f"expected {expected_src_sha}, got {actual_src_sha}"
        )

    dm = _import_dist_m4ri(repo_root)
    stim = dm._get_stim()
    stim_ver = getattr(stim, "__version__", "unknown")
    req_stim = recipe.get("stim_version")
    if req_stim and stim_ver != req_stim and not allow_stim_mismatch:
        raise RuntimeError(
            f"Stim version mismatch: manifest pins {req_stim}, installed {stim_ver}. "
            f"Install stim=={req_stim} or pass --allow-stim-mismatch to inspect."
        )

    t0 = time.monotonic()
    circuit = stim.Circuit.from_file(str(src_path))
    if recipe.get("remove_empty_detectors", True):
        circuit, _ = dm.remove_empty_detectors(circuit)
    if recipe.get("rounds") is not None:
        circuit = dm.set_circuit_rounds(circuit, int(recipe["rounds"]))

    noise_spec = recipe.get("add_noise")
    if noise_spec:
        when = noise_spec.get("when", "if_noiseless")
        p_noise = float(noise_spec["p"])
        if when == "always":
            circuit = dm.add_noise(circuit, p=p_noise, force=True)
        elif when == "if_noiseless":
            if not dm.has_noise(circuit):
                circuit = dm.add_noise(circuit, p=p_noise)
        else:
            raise ValueError(f"Unknown add_noise.when={when!r}")

    decomp = bool(recipe.get("decompose_errors", False))
    approx = bool(recipe.get("approximate_disjoint_errors", False))
    dem = circuit.detector_error_model(
        decompose_errors=decomp,
        approximate_disjoint_errors=approx,
    )
    data, meta = canonical_dem_bytes(dem, dm)
    meta["sha256"] = sha256_bytes(data)
    meta["bytes"] = len(data)
    meta["gen_s"] = round(time.monotonic() - t0, 4)
    meta["stim_version"] = stim_ver
    meta["dist_m4ri_py_version"] = getattr(dm, "__version__", "unknown")
    return data, meta


def cache_path_for_sha(cache_dir: Path, sha256: str, ext: str = ".dem") -> Path:
    """Content-addressed storage path under cache_dir."""
    return cache_dir / "sha256" / sha256[:2] / f"{sha256}{ext}"


def safe_cache_name(inst_id: str) -> str:
    """Filesystem-safe name derived from an instance id."""
    return inst_id.replace("/", "__").replace(":", "_")


def write_cache_entry(cache_dir: Path, inst_id: str, data: bytes, sha: str, ext: str = ".dem") -> Path:
    """Write content-addressed cache entry atomically and refresh by_id symlink."""
    dst = cache_path_for_sha(cache_dir, sha, ext=ext)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.is_file():
        tmp = dst.with_name(f".{dst.name}.{os.getpid()}.tmp")
        tmp.write_bytes(data)
        os.replace(tmp, dst)
    else:
        existing_sha = sha256_file(dst)
        if existing_sha != sha:
            raise ValueError(
                f"Corrupted cache file {dst}: expected sha256 {sha}, got {existing_sha}"
            )

    link_dir = cache_dir / "by_id"
    link_dir.mkdir(parents=True, exist_ok=True)
    link_path = link_dir / f"{safe_cache_name(inst_id)}{ext}"
    rel_target = os.path.relpath(dst, link_dir)
    tmp_link = link_dir / f".{link_path.name}.{os.getpid()}.tmp"
    if tmp_link.is_symlink() or tmp_link.exists():
        tmp_link.unlink()
    os.symlink(rel_target, tmp_link)
    os.replace(tmp_link, link_path)
    return dst


def resolve_instance_files(
    inst: Dict[str, Any],
    repo_root: Path = REPO_ROOT,
    cache_dir: Optional[Path] = None,
    verify_sha: bool = True,
) -> Dict[str, Path]:
    """Resolve and optionally sha256-verify every input file of a manifest instance."""
    cdir = cache_dir if cache_dir is not None else default_cache_dir(repo_root)
    resolved: Dict[str, Path] = {}
    for role, spec in inst.get("files", {}).items():
        expected_sha = spec.get("sha256")
        if spec.get("generated"):
            if not expected_sha:
                raise ValueError(
                    f"Instance {inst['id']} file {role!r} is generated but has no sha256 in manifest"
                )
            ext = Path(spec.get("cache_name", ".dem")).suffix or ".dem"
            path = cache_path_for_sha(cdir, expected_sha, ext=ext)
            if not path.is_file():
                raise FileNotFoundError(
                    f"Generated input missing for {inst['id']} ({role}): {path}. "
                    f"Run benchmark/manifest/make_inputs.py first."
                )
        else:
            path = repo_root / spec["path"]
            if not path.is_file():
                raise FileNotFoundError(f"Tracked input missing for {inst['id']} ({role}): {path}")
        if verify_sha and expected_sha:
            actual_sha = sha256_file(path)
            if actual_sha != expected_sha:
                raise ValueError(
                    f"sha256 mismatch for {inst['id']} ({role}: {path}): "
                    f"expected {expected_sha}, got {actual_sha}"
                )
        resolved[role] = path
    return resolved


def select_instances(
    manifest: Dict[str, Any],
    tiers: Optional[Sequence[str]] = None,
    ids: Optional[Sequence[str]] = None,
    only_small: bool = False,
) -> List[Dict[str, Any]]:
    """Filter manifest instances by tier, id, and small/large DEM flag."""
    tier_set = set(tiers) if tiers else None
    id_set = set(ids) if ids else None
    out: List[Dict[str, Any]] = []
    for inst in manifest["instances"]:
        if tier_set is not None and inst["tier"] not in tier_set:
            continue
        if id_set is not None and inst["id"] not in id_set:
            continue
        if only_small and inst.get("kind") == "dem" and (inst.get("recipe") or {}).get("large"):
            continue
        out.append(inst)
    return out


def verify_manifest_inputs(
    manifest: Dict[str, Any],
    instances: Iterable[Dict[str, Any]],
    repo_root: Path = REPO_ROOT,
    cache_dir: Optional[Path] = None,
    warm_page_cache: bool = False,
) -> Dict[str, Any]:
    """Verify sha256 of all files of the selected instances (and optionally warm OS page cache)."""
    cdir = cache_dir if cache_dir is not None else default_cache_dir(repo_root)
    inst_list = list(instances)
    checked_files: Dict[Path, int] = {}
    errors: List[str] = []
    for inst in inst_list:
        try:
            files = resolve_instance_files(inst, repo_root=repo_root, cache_dir=cdir, verify_sha=True)
            for p in files.values():
                if p not in checked_files:
                    checked_files[p] = p.stat().st_size
                    if warm_page_cache:
                        # sha256_file already read the full file into the OS page cache.
                        pass
        except Exception as exc:
            errors.append(f"{inst['id']}: {exc}")
    if errors:
        raise RuntimeError(
            "Input verification failed:\n  " + "\n  ".join(errors)
        )
    return {
        "checked": len(inst_list),
        "files": len(checked_files),
        "bytes": sum(checked_files.values()),
        "sha256_ok": True,
        "page_cache_warmed": bool(warm_page_cache),
    }


def _col_multiset(path_or_text: Any, dm: Any) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    cols = dm._dem_error_columns(path_or_text)
    return sorted((tuple(sorted(s)), tuple(sorted(o))) for s, o in cols if s or o)


def _line_stats(path: Path) -> Dict[str, int]:
    n_err = n_det = n_obs = n_other = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b"error("):
                n_err += 1
            elif line.startswith(b"detector"):
                n_det += 1
            elif line.startswith(b"logical_observable"):
                n_obs += 1
            elif line.strip():
                n_other += 1
    return {"error": n_err, "detector": n_det, "logical_observable": n_obs, "other": n_other}


# Known existing DEM files in scratch/ that overlap with certify_v1.json instances.
OVERLAP_MAP: Dict[str, List[Tuple[str, str]]] = {
    "dem/bb_72_12_6_si1000_r6_X.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/bb_72_12_6_si1000_r6_X_full.dem"),
    ],
    "dem/mitten_150_30_8_C5xS3_custom_r5_Z.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/mitten_150_30_8_C5xS3_custom_r5_Z_full.dem"),
    ],
    "dem/bb_90_8_10_si1000_r10_X.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/bb_90_8_10_si1000_r10_X_full.dem"),
    ],
    "dem/torus_n60_k8_d14_X.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/torus_n60_k8_d14_X_full.dem"),
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/torus_n60_k8_d14_X_cur_full.dem"),
    ],
    "dem/mitten_270_54_10_C9xS3_hef_r4_Z.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/mitten_270_54_10_C9xS3_hef_r4_Z_full.dem"),
    ],
    "dem/shyps_225_16_8_r9_Z.full": [
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/shyps_225_16_8_r9_Z_full.dem"),
    ],
    "dem/torus_n200_k8_d17_X.full": [
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/torus_n200_k8_d17_X_full.dem"),
    ],
    "dem/gross_144_12_12_si1000_r12_X.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/gross_144_12_12_si1000_r12_X_full.dem"),
    ],
    "dem/gross_144_12_12_uniform_r12_X.full": [
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/gross_144_12_12_uniform_r12_X_full.dem"),
    ],
    "dem/torus_n200_k16_d16_X.full": [
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/torus_n200_k16_d16_X_full.dem"),
    ],
    "dem/yarn_mitten_150_30_10_hookfree_r10_X.full": [
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/yarn_mitten_150_30_10_hookfree_r10_X_full.dem"),
    ],
    "dem/yarn_mitten_150_30_10_hookfree_r10_Z.full": [
        ("bench_0.15.1_p0", "scratch/bench_0.15.1_p0/work/dems/yarn_mitten_150_30_10_hookfree_r10_Z_full.dem"),
    ],
    "dem/mitten_270_54_10_C9xS3_custom_r4_Z.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/mitten_270_54_10_C9xS3_custom_r4_Z_full.dem"),
    ],
    "dem/double_gross_288_12_18_si1000_r18_X.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/double_gross_288_12_18_si1000_r18_X_full.dem"),
    ],
    "dem/double_gross_288_si1000_X.T1": [
        ("bench_0.17", "scratch/bench_0.17/dems/dgross_288_T1.dem"),
    ],
    "dem/bb_360_12_24_si1000_r24_X.full": [
        ("bench_0.14.0", "scratch/bench_0.14.0/data/bb_360_12_24_si1000_r24_X_full.dem"),
    ],
}


def check_overlap(
    manifest: Dict[str, Any],
    cache_dir: Path,
    scratch_root: Path,
    repo_root: Path = REPO_ROOT,
) -> List[Dict[str, Any]]:
    """Compare cached canonical DEMs against existing DEMs in scratch/bench_0.14.0, 0.15.1_p0, 0.17."""
    dm = _import_dist_m4ri(repo_root)
    report: List[Dict[str, Any]] = []
    by_id = {inst["id"]: inst for inst in manifest["instances"]}
    for inst_id, refs in OVERLAP_MAP.items():
        inst = by_id.get(inst_id)
        if not inst:
            continue
        dem_spec = inst["files"].get("dem", {})
        expected_sha = dem_spec.get("sha256")
        cached = cache_path_for_sha(cache_dir, expected_sha, ".dem") if expected_sha else None
        is_cached = cached is not None and cached.is_file()
        for src_label, rel_path in refs:
            ref_path = scratch_root / rel_path
            entry: Dict[str, Any] = {
                "id": inst_id,
                "source": src_label,
                "ref_path": str(ref_path),
                "ref_exists": ref_path.is_file(),
                "cached": is_cached,
                "expected_sha256": expected_sha,
            }
            if not ref_path.is_file():
                entry["status"] = "ref_missing"
                report.append(entry)
                continue
            ref_sha = sha256_file(ref_path)
            entry["ref_sha256"] = ref_sha
            if not is_cached:
                entry["status"] = "deferred_not_cached"
                entry["manifest_sha_matches_ref"] = (ref_sha == expected_sha)
                report.append(entry)
                continue
            assert cached is not None
            gen_sha = sha256_file(cached)
            entry["gen_sha256"] = gen_sha
            if gen_sha == ref_sha:
                entry["status"] = "exact_sha256_match"
                entry["diff_summary"] = None
            else:
                gen_lines = _line_stats(cached)
                ref_lines = _line_stats(ref_path)
                # Compare error lines verbatim and as GF(2) column multisets
                gen_err = [ln for ln in cached.read_text().splitlines() if ln.startswith("error(")]
                ref_err = [ln for ln in ref_path.read_text().splitlines() if ln.startswith("error(")]
                err_lines_equal = (gen_err == ref_err)
                cols_equal = err_lines_equal or (_col_multiset(str(cached), dm) == _col_multiset(str(ref_path), dm))
                entry["status"] = "format_diff_same_errors" if err_lines_equal else (
                    "order_or_prob_diff_same_cols" if cols_equal else "matrix_mismatch"
                )
                entry["gen_lines"] = gen_lines
                entry["ref_lines"] = ref_lines
                entry["error_lines_identical"] = err_lines_equal
                entry["columns_identical"] = cols_equal
                diff_parts = []
                if ref_lines["detector"] != gen_lines["detector"]:
                    diff_parts.append(
                        f"ref has {ref_lines['detector']} detector() coordinate lines vs "
                        f"{gen_lines['detector']} in canonical DEM"
                    )
                if ref_path.read_bytes().endswith(b"\n") != cached.read_bytes().endswith(b"\n"):
                    diff_parts.append("trailing newline differs")
                if not err_lines_equal:
                    diff_parts.append(f"error lines differ (columns_identical={cols_equal})")
                entry["diff_summary"] = "; ".join(diff_parts) or "byte difference"
            report.append(entry)
    return report


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate and sha256-verify benchmark inputs from a manifest."
    )
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Path to manifest JSON.")
    p.add_argument("--repo-root", type=Path, default=REPO_ROOT, help="Repository root containing benchmark/.")
    p.add_argument("--cache-dir", type=Path, default=None, help="Content-addressed cache directory.")
    p.add_argument("--tier", type=str, default=None, help="Comma-separated tiers to include (e.g. T0,T1).")
    p.add_argument("--id", action="append", default=None, help="Specific instance id(s) to process.")
    p.add_argument("--only-small", action="store_true", help="Skip large DEMs (recipe.large=true).")
    p.add_argument("--verify-only", action="store_true", help="Verify sha256s without generating missing DEMs.")
    p.add_argument("--force", action="store_true", help="Regenerate DEMs even if already in cache.")
    p.add_argument("--allow-stim-mismatch", action="store_true", help="Allow a different stim version (inspection).")
    p.add_argument("--check-overlap", action="store_true", help="Compare cached DEMs against scratch/ reference DEMs.")
    p.add_argument(
        "--scratch-root",
        type=Path,
        default=REPO_ROOT,
        help="Root directory holding scratch/bench_0.14.0 and scratch/bench_0.15.1_p0 (for --check-overlap).",
    )
    p.add_argument("--json-out", type=Path, default=None, help="Write machine-readable summary JSON to this file.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    repo_root = args.repo_root.resolve()
    cache_dir = (args.cache_dir if args.cache_dir is not None else default_cache_dir(repo_root)).resolve()
    manifest, manifest_sha = load_manifest(args.manifest.resolve())
    tiers = [t.strip() for t in args.tier.split(",") if t.strip()] if args.tier else None
    instances = select_instances(manifest, tiers=tiers, ids=args.id, only_small=args.only_small)
    skipped_large = [
        inst["id"] for inst in select_instances(manifest, tiers=tiers, ids=args.id, only_small=False)
        if args.only_small and inst.get("kind") == "dem" and (inst.get("recipe") or {}).get("large")
    ]

    print(
        f"Manifest: {manifest['id']} v{manifest['version']} (sha256={manifest_sha[:16]}), "
        f"selected={len(instances)}, skipped_large={len(skipped_large)}, cache={cache_dir}",
        flush=True,
    )

    generated: List[Dict[str, Any]] = []
    verified_cached = 0
    verified_tracked = 0
    mismatches: List[str] = []

    for inst in instances:
        inst_id = inst["id"]
        for role, spec in inst.get("files", {}).items():
            expected_sha = spec.get("sha256")
            if spec.get("generated"):
                ext = Path(spec.get("cache_name", ".dem")).suffix or ".dem"
                cached = cache_path_for_sha(cache_dir, expected_sha, ext=ext) if expected_sha else None
                if not args.force and cached is not None and cached.is_file():
                    actual = sha256_file(cached)
                    if actual != expected_sha:
                        mismatches.append(
                            f"{inst_id} ({role}): cached sha256 mismatch expected={expected_sha} got={actual}"
                        )
                    else:
                        verified_cached += 1
                    continue
                if args.verify_only:
                    mismatches.append(f"{inst_id} ({role}): generated file missing from cache ({cached})")
                    continue
                recipe = inst.get("recipe")
                if not recipe:
                    mismatches.append(f"{inst_id} ({role}): missing recipe for generated file")
                    continue
                data, meta = generate_canonical_dem(
                    recipe, repo_root=repo_root, allow_stim_mismatch=args.allow_stim_mismatch
                )
                actual = meta["sha256"]
                if expected_sha and actual != expected_sha:
                    mismatches.append(
                        f"{inst_id} ({role}): generated sha256 mismatch "
                        f"expected={expected_sha} got={actual} (stim={meta['stim_version']})"
                    )
                    continue
                dst = write_cache_entry(cache_dir, inst_id, data, actual, ext=ext)
                generated.append({"id": inst_id, "path": str(dst), **meta})
                print(
                    f"  generated {inst_id}: {meta['detectors']}x{meta['errors']} "
                    f"(k={meta['observables']}, {meta['bytes']} B, {meta['gen_s']:.2f}s) "
                    f"sha256={actual[:16]} OK",
                    flush=True,
                )
            else:
                path = repo_root / spec["path"]
                if not path.is_file():
                    mismatches.append(f"{inst_id} ({role}): missing tracked file {path}")
                    continue
                actual = sha256_file(path)
                if expected_sha and actual != expected_sha:
                    mismatches.append(
                        f"{inst_id} ({role}): tracked sha256 mismatch {spec['path']} "
                        f"expected={expected_sha} got={actual}"
                    )
                    continue
                verified_tracked += 1

    if skipped_large:
        print("Deferred large DEMs (--only-small):", flush=True)
        for sid in skipped_large:
            print(f"  - {sid}", flush=True)

    overlap_report: List[Dict[str, Any]] = []
    if args.check_overlap:
        overlap_report = check_overlap(manifest, cache_dir, args.scratch_root.resolve(), repo_root=repo_root)
        print("Overlap check against existing scratch/ DEMs:", flush=True)
        for row in overlap_report:
            extra = f" ({row['diff_summary']})" if row.get("diff_summary") else ""
            if row["status"] == "deferred_not_cached":
                extra = f" (manifest_sha_matches_ref={row['manifest_sha_matches_ref']})"
            print(f"  {row['id']:45s} [{row['source']:13s}] -> {row['status']}{extra}", flush=True)

    summary = {
        "manifest": {"id": manifest["id"], "version": manifest["version"], "sha256": manifest_sha},
        "cache_dir": str(cache_dir),
        "selected": len(instances),
        "verified_tracked": verified_tracked,
        "verified_cached": verified_cached,
        "generated": generated,
        "skipped_large": skipped_large,
        "mismatches": mismatches,
        "overlap": overlap_report,
    }
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(summary, indent=2) + "\n")

    if mismatches:
        print("ERROR: sha256 / input verification failed:", file=sys.stderr)
        for msg in mismatches:
            print(f"  - {msg}", file=sys.stderr)
        return 1

    print(
        f"OK: verified {verified_tracked} tracked files, {verified_cached} cached DEMs, "
        f"generated {len(generated)} DEMs, 0 mismatches.",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
