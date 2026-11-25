---
geometry: "top=0.5in, bottom=0.5in, left=0.25in, right=0.25in"
header-includes:
  - |
    \makeatletter
    \DeclareRobustCommand{\_}{\ifmmode\nfss@text{\textunderscore}\else\textunderscore\allowbreak\fi}
    \makeatother
    \setlength{\emergencystretch}{1.5em}
---

# Benchmark Manifests (`benchmark/manifest/`)

This directory holds versioned instance manifests for the `dist-m4ri` certification benchmark (`tools/bench/`) and
the deterministic input generator `make_inputs.py`. The design follows `scratch/bench_review/REVIEW.md`, Section 3.2.

## Files

- `certify_v1.json`: manifest `certify` v1.0 (51 instances from `benchmark/`, tiers T0-T3).
- `make_inputs.py`: generates canonical DEMs into a content-addressed cache outside git (default
  `scratch/bench_cache/`, or `DM4RI_BENCH_CACHE`) and verifies every sha256. Refuses to proceed on any mismatch or
  unpinned Stim version.

## Manifest schema (`dm4ri-bench-manifest/1`)

Top-level object:

- `schema`: `"dm4ri-bench-manifest/1"`.
- `id`, `version`, `created`, `description`: manifest identity (tiers are frozen within a version).
- `tiers`: map `"T0".."T3"` to `{range_s: [lo, hi|null], budget_s, seeds: [int], threads: [int]}`:
  - `T0` (< 1 s): budget 60 s, seeds 1..5, 128 threads;
  - `T1` (1-60 s): budget 300 s, seeds 1..5, 128 and 64 threads;
  - `T2` (1-30 min): budget 1800 s, seeds 1..3, 128 and 64 threads;
  - `T3` (> 30 min / open): budget 3600 s, seed 1, 128 threads.
- `canary_instance`: fixed T1 instance run periodically by the harness to detect performance drift.
- `scaling_instances`: Track S thread-scaling inputs (`method=2` fixed work and `method=3`).
- `instances`: list of instance objects.

Instance object:

- `id`: unique identifier:
  - `css/<name>:Z` (`finH=Hx finG=Hz`, giving $d_Z$) or `css/<name>:X` (`finH=Hz finG=Hx`, giving $d_X$);
  - `classical/<name>` (`finH=H classical=1`);
  - `hl/<name>` (`finH=H finL=L`, optional for future manifests);
  - `dem/<name>.full` or `dem/<name>.T1` (generated canonical DEM) or `dem/<name>` (tracked `.dem`);
  - `circ/<name>` (tracked `.stim` circuit; engines that do not read `.stim` convert it via `Adapter.prepare()`).
- `kind`: `"css_sector"`, `"classical"`, `"hl"`, `"dem"`, or `"circuit"`.
- `tier`: `"T0"`, `"T1"`, `"T2"`, or `"T3"`.
- `sector`: `"Z"`, `"X"`, or `null`.
- `d_ref`: known exact distance (`int` or `null`), and `bracket_ref`: `[lo, hi]` or `null` when exact $d$ is open.
- `ref_source`: list of references (sections of `BENCHMARK.md` / `HISTORY.md`, engine, date).
- `t_best_s`, `t_best_note`: baseline calibration time (seconds, or `null`) and provenance note.
- `compare3_id`: optional external benchmark cross-reference id, or `null`.
- `tags`: `{n, m, k, family, ...}` (for DEMs/circuits: `m` detectors, `n` errors in the canonical DEM, `k`
  observables, plus `code_n`, `code_k`, `rounds`, `basis`, `noise`).
- `files`: map of role (`H`, `G`, `L`, `dem`, `circuit`) to either:
  - tracked file: `{path: "benchmark/...", sha256: "..."}`, or
  - generated file: `{generated: true, cache_name: "...", sha256: "..."}`, stored at
    `<cache_dir>/sha256/<sha[:2]>/<sha>.dem` (with a convenience symlink in `<cache_dir>/by_id/`).
- `recipe`: `null` for tracked matrices/DEMs, or `canonical_dem_v1`:
  - `source_stim`: `{path, sha256}` of the tracked `.stim` file;
  - `stim_version`: pinned Stim version (`"1.16.0"`);
  - `dist_m4ri_py_version`: `dist_m4ri.py` version whose `add_noise` / `_dem_error_columns` was used (`"0.17.3"`);
  - `remove_empty_detectors`: `true`;
  - `rounds`: `null` (keep file's repeat count) or integer (e.g. `0` for `.T1`);
  - `add_noise`: `null` or `{p: 0.001, when: "if_noiseless"}` (for noiseless torus templates);
  - `decompose_errors`: `false`, `approximate_disjoint_errors`: `false`;
  - `canonical_format`: `"run_bench_0.14.0"` — one `error(<p>) D... L...` line per non-empty mechanism (sorted
    targets, `^` components summed mod 2, no `detector(...)` coordinate lines), plus `detector D<m-1>` /
    `logical_observable L<k-1>` only if the last index is absent from all errors;
  - `large`: `true` for multi-megabyte DEMs that `--only-small` skips.

## Usage

```bash
# Generate and sha256-verify all inputs (uses scratch/bench_cache/ by default)
python3 benchmark/manifest/make_inputs.py

# Generate only small DEMs (defers gross, double_gross, bb_360 full DEMs)
python3 benchmark/manifest/make_inputs.py --only-small

# Verify sha256 of cached + tracked inputs for a tier without generating anything
python3 benchmark/manifest/make_inputs.py --verify-only --tier T0

# Cross-check cached canonical DEMs against existing scratch/ DEMs
python3 benchmark/manifest/make_inputs.py --only-small --check-overlap
```
