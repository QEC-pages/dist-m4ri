# Record Schemas of the Certification Benchmark (`tools/bench/`)

This file defines the JSON written by `harness.py` and read by `verify.py` and `analyze.py`: the run record
`dm4ri-bench-run/1`, the environment `dm4ri-bench-env/1`, and the plan `dm4ri-bench-plan/1`. The instance manifest
`dm4ri-bench-manifest/1` is described in [`benchmark/manifest/README.md`](../../benchmark/manifest/README.md). The
design is that of `scratch/bench_review/REVIEW.md`, Sections 3.3-3.5. The field names follow qec-cert's `compare3`
(`qec-cert/bench/compare3/results.json`, schema `qec-cert-compare3-results/v1`) where it has them: `lb`, `ub`,
`exact`, `wall_sec`, `cap_sec`, `threads`, `lock_wait_sec`, `witness`, `witness_weight`, `witness_verified_numpy`,
`witness_check_detail`, and the `tools` block (here: `engine`).

## 1. Session directory

```
<session>/
  session.json     # session id, harness version, command lines, options, input check (dm4ri-bench-session/1)
  env.json         # the environment (Section 4)
  plan.json        # the planned runs in their run order (Section 5)
  results.jsonl    # one run record per line (Section 2), appended atomically; the LAST record of a run_id wins
  runs/<run_id>/   # run.json (the same record), out.txt (stdout), err.txt (stderr)
  harness.log      # one line per event
```

## 2. The run record `dm4ri-bench-run/1`

One JSON object per run. Times are in seconds, sizes `[m, n]` are rows x columns, and every bound is an integer. A
field that does not apply, or is unknown, is `null`; `harness.py` writes every field listed here.

**Top level**
- `schema`: `"dm4ri-bench-run/1"`.
- `run_id`: sha1 (hex) of the canonical JSON (sorted keys, no spaces) of `key`.
- `key`: `manifest` (`"certify/1.0"`), `instance`, `engine`, `config`, `league`, `threads`, `budget_s`, `seed`,
  `canary` (int or null), `harness_major` (1).
- `session`: e.g. `"2026-11-25T01:00:00Z-betamach"`.
- `harness`: `version` (e.g. `"1.0.0"`).
- `manifest`: `id`, `version`, `sha256` (of the manifest file).
- `instance`: `id`, `kind`, `tier`, `d_ref`, `bracket_ref` (`[lo, hi]` or null), `compare3_id`.
- `engine`: the description of the engine (Section 3).
- `options`: `league`, `threads`, `budget_s`, `seed`, `argv_extra` (list of str), `hints` (object), `canary` (int
  or null), `kill_after_s`, `cpus` (the CPU list of `taskset`). The leagues: `default` (only `threads`, the budget,
  and the seed are set: Track C), `scaling` (Track S: the binary's thread scaling, configurations `default` and
  `cc_fixed`), and `tuned` (any options and hints).
- `placement`: `taskset` (CPU list or null), `nice`, `exclusive` (bool: the lock was held and the load gate
  passed), `lock` (path or null), `lock_wait_sec`, `load_gate` {`busy_cores`, `max_busy`, `waited_s`, `passed`}.
- `result`, `verify`, `env`: below.
- `raw`: `dir`, `stdout`, `stderr` (paths relative to the session directory), `argv` (list of str).

**`result`**
- `status`: the final status, `certified`, `bracket`, `timeout`, `error`, or `wrong` (after `verify`).
- `status_engine`: the status from the engine's output alone (never `wrong`).
- `lb`, `ub`: the bounds the engine reported (a killed run: the last ones it printed); `ub` null: none found.
- `exact`: `lb == ub > 0`, as the engine reported.
- `wall_sec`: the primary time: the wall time of the process from fork to exit, measured by the harness. It includes
  the start-up, the input stage, and (Python) the imports.
- `cap_sec`: the budget T given to the engine (= `options.budget_s`).
- `prep_s`: the time to convert the input into the engine's format (0 for the native formats).
- `cpu_user_s`, `cpu_sys_s`, `max_rss_kb`: CPU time and maximum resident set size (kB) of the process and its reaped
  children, from `os.wait4`.
- `rc` (exit code; null if a signal ended it), `signal`, `killed` (the harness killed it at `options.kill_after_s` =
  T + max(30, 0.05 T): SIGTERM, then SIGKILL 10 s later).
- `t_internal`: the engine's own times: `t_stop` (binary: the search time, from the end of the input stage),
  `t_input` (binary: process start to search start, if reported), `t_import`, `t_call` (Python).
- `stop_reason` (the stop line, as text) and `stop_code`, a stable token: `bounds`, `cc_found`, `cc_done`,
  `cc_slow`, `timeout`, `rw_converged`, `rw_steps`, `wmin`, `maxc`, `collected`, `dstop`, `k0`, `none`, `killed`,
  `other`.
- `t_ub`, `t_lb`: the time (same clock as `t_stop`) at which the final ub / lb was reached, if the engine stamps it.
- `rw_steps`: RW information sets (binary).
- `size`: `full` `[m, n]` (as read), `searched` `[m, n]` (after merging and dropping).
- `cc_rounds`: per CC round {`w`, `thread_s`, `wall_s`, `nodes`, `complete`} (empty if unknown).
- `threads_used`: `requested`, `rw`, `cc`.
- `rw_opt`: the resolved RW stage mask; `k0`: the code has no logical operators (k = 0).
- `witness`: null, or {`weight`, `support` (0-based columns of the input as given to the engine), `base` (always 0),
  `source` (`rw`, `cc`, `finc`, `cert`, ...), `complete` (false if the engine printed it truncated or garbled)}.
- `events`: the bounds over time, a list of `[t, lb, ub]` (null where unknown).
- `parse`: how the result was read: `json` (the binary's `# RESULT` line), `stderr` (its text), `worker` (Python).
- `engine_info`: engine-specific extras (e.g. `engine_auto` of `dist_m4ri.py`).
- `error`: null, or {`kind`: `engine` or `harness`, `msg`}.

**Status rules** (the adapter, then `verify.py`)
- `certified`: the engine reported `lb == ub > 0` before it ended, and no check failed.
- `bracket`: the engine ended by itself before the budget with `lb < ub` (e.g. the binary's `cc_slow` stop).
- `timeout`: the engine stopped at its own time limit, or the harness killed it; its bracket is the last one reported.
- `error`: a nonzero exit code, a crash, or an unreadable output (`error.kind = engine`), or a failure of the harness
  itself (`error.kind = harness`; `--retry` runs such runs again).
- `wrong`: a check of `verify.py` failed (Section 2.1). It counts as unsolved and is listed in Table C.

A run is **solved** if its final status is `certified`. PAR-2 scores a solved run by its `wall_sec`, and any other
run by 2 T (T = `options.budget_s`).

**`verify`** (written by `verify.py`, outside the timed region)
- `witness_weight`: the weight of the witness checked (null: none).
- `witness_verified_numpy`: the GF(2) check passed (null: no witness, or an incomplete one).
- `witness_check_detail`: as in compare3: `ok`, `weight`, `H_syndrome_weight`, `L_syndrome_weight`,
  `independent_of_G`, `reason`.
- `ub_verified`: this run has a verified witness of weight `ub`.
- `lb_audit`: `ok`, `violation`, or `n/a` (this run's lb against `d_ref`/`bracket_ref` and its own witness).
- `matches_ref`: a certified value equals `d_ref` (or lies in `bracket_ref`); null if not certified or no reference.
- `violations`: the text of each failed check (not empty only for the status `wrong`).
- `unverified_reason`: why the upper bound could not be verified (no witness printed, a truncated support, ...).

### 2.1 Checks that make a run `wrong`
1. The witness fails the GF(2) check: `H e != 0`, or `e` is trivial (`L e = 0`, or `e` in the row space of `G`).
2. A verified witness is lighter than the `lb` of the same run, or than `d_ref` (then the reference is wrong too).
3. `lb > d_ref`, or `lb > bracket_ref[1]`.
4. A certified value differs from `d_ref`, or lies outside `bracket_ref`.

The cross-run audit (`verify.py --audit`, and Table C) adds: the `lb` of any run above the weight of any verified
witness of the same instance (from any engine or run).

**`env`** (per run)
- `env_id`: sha1 of the canonical JSON of the session's `env.json`.
- `busy_cpu_s`: the busy CPU time of the machine over the run (`/proc/stat`: user, nice, system, irq, softirq, steal).
- `child_cpu_s`: `cpu_user_s + cpu_sys_s`; `harness_cpu_s`: the CPU time of the harness process over the run.
- `ext_cpu_s`: `busy_cpu_s - child_cpu_s - harness_cpu_s`, the CPU time of other processes during the run.
- `steal_s`: the steal time of the machine over the run.
- `contaminated`: `ext_cpu_s > 0.02 * threads * wall_sec`.
- `started_at`, `ended_at`: UTC, ISO 8601.

## 3. The engine block (`Adapter.describe()`)

- `name`: `dist_m4ri` (the binary), `dist_m4ri.py` (Python), `qec_cert`, `tropidist`.
- `interface`: `binary`, `python`, `library`, `cli`.
- `config`: `default` (the binary, `method=3`), `cc_fixed` (the binary, `method=2 wmax=d_ref-1`: fixed work, Track S
  only), `off`, `auto`, `cert` (Python `engine=`).
- `version`, `commit`, `dirty`, `sha256`: the engine version, git commit, uncommitted changes, and the sha256 of the
  binary, library, or module.
- `deps`: e.g. the versions of `m4ri`, `stim`, `numpy`, `scipy`, `qec_cert`.
- `features`: e.g. `{"json": true}`: the binary's machine-readable output (`json=1 stamp=1`) was used.

An **engine entry** of the analysis is `name:config@p<threads>`, e.g. `dist_m4ri:default@p128` (the baseline) or
`dist_m4ri.py:auto@p128`. Only records of the league `default` that are not canaries (`options.canary` null) enter
Tables A-B and Figures 1-4. Figure 5 uses the league `scaling` (and the `default` records of the binary at p = 64 and
128); the canaries are shown as a drift check.

## 4. `env.json` (`dm4ri-bench-env/1`)

- `schema`, `created` (UTC), `host`, `kernel`, `os`, `vm`, `sched` {`autogroup`}.
- `cpu`: `model`, `sockets`, `cores_per_socket`, `threads_per_core`, `cpus`, `numa`, `caches`, `flags_sha1`,
  `hypervisor`, `lscpu` (the raw text).
- `memory`: `total_kb`.
- `compiler`: `cc` (version), `flags` (from `src/makefile`); `m4ri` (version).
- `python`: `version`, `executable`, `stim`, `numpy`, `scipy`, `matplotlib`.
- `qec_cert`: `path`, `version`, `commit`, `dirty`, `lib`, `lib_sha256`.
- `dist_m4ri`: `repo`, `version`, `commit`, `dirty`, `binary`, `sha256`, `py`, `py_sha256`, `py_version`.
- `harness`: `version`, `files` (the sha256 of each script).
- `manifest`: `id`, `version`, `sha256`.
- `inputs`: `checked`, `files`, `bytes`, `sha256_ok`, `page_cache_warmed`.

## 5. `plan.json` (`dm4ri-bench-plan/1`)

`schema`, `created`, `order_seed`, `manifest`, `options` (the tiers, budgets, seeds, threads, and engines), and `runs`:
a list in run order of {`run_id`, `key`, `tier`, `instance`, `engine`, `config`, `threads`, `budget_s`, `seed`}.
Within a tier the runs are a random permutation (from `order_seed`); the tiers run in the order T0, T1, T2, T3. The
canary runs are not in the plan: the harness inserts one when the last canary is older than `--canary-every` (default
1800 s).
