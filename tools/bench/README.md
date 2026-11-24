# Certification Benchmark Harness (`tools/bench/`)

Engine-ready certification benchmark for `dist-m4ri` (design: `scratch/bench_review/REVIEW.md`, Sections 3-4;
record schemas: [`SCHEMA.md`](SCHEMA.md); manifest schema:
[`benchmark/manifest/README.md`](../../benchmark/manifest/README.md)).

## Layout

- `benchmark/manifest/certify_v1.json`: frozen v1.0 manifest (51 instances in tiers T0-T3, sha256s, DEM recipes,
  `d_ref`/`bracket_ref`, 22 `compare3` links).
- `benchmark/manifest/make_inputs.py`: deterministic canonical DEM generation into a content-addressed cache
  outside git (default `scratch/bench_cache/`, or `DM4RI_BENCH_CACHE`) and sha256 verification.
- `tools/bench/harness.py`: randomized per-tier plan, exclusive runs under `flock /tmp/qec-cert-bench.lock`
  (shared with `qec-cert`), load gate, `os.wait4` CPU/RSS + `/proc/stat` contamination recording, periodic canary
  runs, `env.json`, and restartable atomic JSONL logging.
- `tools/bench/adapters/`: uniform `Adapter` interface for `dist_m4ri` binary (`default` and `cc_fixed`),
  `dist_m4ri.py` (`off`, `auto`, `cert`), plus stubs for `qec_cert` and `tropidist` (`supports()` False until
  configured).
- `tools/bench/verify.py`: independent numpy GF(2) witness check and lower-bound cross-audit (marks violations
  `"wrong"`).
- `tools/bench/analyze.py`: pure function of manifest + records + `env.json` generating `REPORT.md`,
  `summary.json`, Tables A-D, and Figures 1-5.

## 1. Generate and verify inputs

```bash
# Generate all canonical DEMs into scratch/bench_cache/ and verify every sha256
python3 benchmark/manifest/make_inputs.py

# Or generate only small DEMs (defers gross, double_gross, bb_360 full DEMs)
python3 benchmark/manifest/make_inputs.py --only-small

# Verify existing cache + tracked files without generating anything
python3 benchmark/manifest/make_inputs.py --verify-only --tier T0,T1
```

## 2. Plan and run a certification session

```bash
# Inspect the full Track C plan without running anything
python3 tools/bench/harness.py --session scratch/bench_sessions/v1_run1 --plan-only

# Run Track C (binary on p=128 + p=64 on T1/T2; budgets 60/300/1800/3600 s; seeds 5/5/3/1)
python3 tools/bench/harness.py --session scratch/bench_sessions/v1_run1 \
    --engines dist_m4ri:default,dist_m4ri.py:off,dist_m4ri.py:auto,dist_m4ri.py:cert

# Run Track S (thread scaling of method=2 fixed work and method=3)
python3 tools/bench/harness.py --session scratch/bench_sessions/v1_scaling \
    --track scaling --engines dist_m4ri:default,dist_m4ri:cc_fixed
```

Pinned smoke test (no lock, no load gate, T0 only, 2 threads):

```bash
taskset -c 126-127 nice -n 19 python3 tools/bench/harness.py \
    --session /tmp/dm4ri_t0_smoke --tiers T0 --threads 2 --seeds 1 --budget 30 \
    --taskset 126-127 --nice 19 --no-lock --no-load-gate --canary-every 0
```

## 3. Restart and retry

Re-running the exact same `--session <dir>` command skips every `run_id` that already has a final record in
`results.jsonl`. Add `--retry` to re-run prior harness-level errors (`error.kind == "harness"`), or
`--retry-all-errors` to also retry engine errors.

## 4. Audit and analyze

```bash
# Cross-run lower-bound and witness audit across one or more sessions
python3 tools/bench/verify.py --manifest benchmark/manifest/certify_v1.json \
    --session scratch/bench_sessions/v1_run1 --audit

# Generate Tables A-D, Figures 1-5, summary.json, and REPORT.md
python3 tools/bench/analyze.py --manifest benchmark/manifest/certify_v1.json \
    --session scratch/bench_sessions/v1_run1 --out scratch/bench_sessions/v1_report
```
