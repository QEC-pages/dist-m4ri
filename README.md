# dist-m4ri - Distance of a Classical or Quantum CSS Code

## Overview

`dist-m4ri` is a high-performance multithreaded C program and Python library for computing and bracketing the minimum
distance of binary classical linear codes, quantum CSS codes, and Stim Detector Error Models (DEMs).

The program implements three main methods (with `method=3` as the default):

- **Method 1 (`method=1`) - Random Window (RW) Algorithm**: Multithreaded random information set search to find
  low-weight non-trivial codewords and establish an **upper distance bound** $d_{\max}$.
- **Method 2 (`method=2`) - Connected Cluster (CC) Algorithm**: Multithreaded exhaustive depth-first cluster
  enumeration to compute **exact distance** or establish a certified **lower distance bound** $d_{\min}$.
- **Method 3 (`method=3`, default) - Bracketing Mode (Artillery Fork / Вилка)**: Concurrently runs CC and RW on multiple
  threads, dynamically balancing CPU cores between CC and RW based on current bounds $[d_{\min}, d_{\max}]$, distance
  estimate (`dexp`/`dest`), remaining RW steps, timeout, and the measured scaling characteristics of CC.

For a classical binary linear code, only the parity-check matrix $H$ is needed.

For a quantum CSS code, matrix $H_X$ and either $H_Z$ (as `finG`) or logical operators $L_X$ (as `finL`) are needed.

Alternatively, a detector error model file from `stim` can be specified using `fdem=[str]`.

All matrices with entries in $\text{GF}(2)$ have $n$ columns and obey the orthogonality conditions:
$$H_X H_Z^T = 0,\quad H_X L_Z^T = 0,\quad L_X H_Z^T = 0,\quad L_X L_Z^T = I.$$

---

## Output Format & Stream Separation

`dist-m4ri` strictly separates machine-parseable results from informational progress logs:

- **`stdout`**: Outputs three space-separated integers:
  ```text
  dmin dmax rw_steps
  ```
  - `dmin - 1` is the maximum cluster size analyzed without success by CC (`dmin = dmax` if CC found a minimum-weight
    codeword).
  - `dmax` is the weight of the smallest non-trivial codeword found by RW (`0` if none found).
  - `rw_steps` is the number of completed RW steps across all threads (`0` if CC found a minimum-weight codeword, or
    if RW did not run in `method=2`).
  - When `dmin = dmax = d`, the exact code distance is confirmed.

> **Note on Compatibility**: This 3-number output format (`dmin dmax rw_steps`) is specific to the multithreaded
> `dist_m4ri` and is incompatible with the legacy single-threaded `dist_m4ri_old` (which returned a single integer `d`
> or `-w`).

- **`stderr`**: Receives all status banners, thread balancing reports, CC round timings, RW discovery logs, warnings,
  and confinement profiles.

---

## How the Methods Work

### 1. Multithreaded RW Algorithm (`method=1`)
Searches for low-weight non-trivial binary codewords $c$ such that $Hc = 0$ and $Lc \neq 0$.
Threads independently generate random column permutations, compute Gaussian elimination to find information sets,
and extract candidate dual-row codewords. When a lighter codeword is discovered, all worker threads atomically update
the global upper bound $d_{\max}$ and prune heavier entries.

Relevant parameters:
- `steps=[int]`: Total number of information sets / RW rounds across all threads (default: 100000).
- `min_hits=[int]`: QDistRnd-style automatic convergence stopping criterion (default: 5; set 0 to disable). Stops RW
  when tracked minimum-weight codewords (up to `cov_cws=100`) have been rediscovered at least `min_hits` times.
- `cov_cws=[int]`: Maximum number of minimum-weight codewords tracked in the hash table for `min_hits` convergence
  (default: 100).
- `kwin=[int]` (alias: `win=[int]`): Localized column permutation window size $W$ (default: 0 for automatic hybrid
  50% uniform / 50% contiguous index window when $n \ge 500$).
- `win_mode=[0|1]`: Window construction mode when `kwin > 0`: `0` for Tanner graph BFS neighborhood, `1` for contiguous
  column index window (default: 0).
- `ksub=[int]`: Subspace dimension sampled from $\ker(H)$ for cache-resident RW elimination (default: 0 for full-matrix
  RW; automatically falls back to `ksub=0` with a warning when $m < \nu = \dim\ker(H)$).
- `refresh=[int]`: RW step interval for adaptive $\ker(H)$ basis refresh via low-weight codeword exchange and
  re-echelonization (default: 5000 when `ksub > 0`, 0 to disable).
- `wmin=[int]`: Minimum distance of interest (stop immediately when a codeword of weight $w \le w_{\min}$ is found).
- `threads=[int]`: Maximum number of POSIX worker threads to run (default: number of CPU cores; subject to automatic
  throttling unless `nothrottle=1` is specified).
- `timeout=[sec]`: Maximum execution time in seconds (default: 60.0).

### 2. Multithreaded CC Algorithm (`method=2`)
Recursively explores connected clusters starting from each column $i \in [0, n-1]$; a cluster started at column $i$
is grown using only columns $j > i$, so that every codeword is found starting from its smallest column. Columns are
distributed dynamically among worker threads via lock-free atomic queues.
If `noscan=0` (default), CC scans weights $w = 1, 2, \dots, w_{\max}$. When `outC` is specified, CC exhausts all
columns for weight $w$ to collect all unique minimum-weight codewords.
The scan ends with the round $w = d$ in which CC finds a codeword (with `outC` and `dW>0`, after the extra rounds up to
$w = d + \text{dW}$). A known upper bound $d_{\max}$ (from `dmax=[int]` or from codewords in `finC`) ends the scan as
soon as $d_{\min} = d_{\max}$, unless codewords are collected (`outC` or `maxC`).

Relevant parameters:
- `wmax=[int]`: Maximum cluster weight to search (optional if `timeout>0` or `dmax>0` is specified; otherwise required
  for CC).
- `smax=[int]`: Maximum syndrome weight to track for confinement profile (default: 0, disabled for faster CC pruning;
  set e.g. `smax=5` to compute confinement).

Expert parameters (each prints a `# WARNING:` to `stderr`; see
[Restricting the CC Search](#restricting-the-cc-search-expert-options) below):
- `noscan=[0|1]`: If set to 1 (`method=2` only), run CC only at weight $w = w_{\max}$, without scanning smaller
  weights. The result is **not** a lower bound unless `dmin=wmax` is also supplied.
- `start=[list]`: Comma-separated list of CC start columns (0-based), e.g. `start=0,48,96`. Clusters grown from a
  listed column are **not** limited to larger column indices (for quasi-cyclic or otherwise symmetric codes).
- `cbeg=[int]`, `cend=[int]`: Range $[c_{\text{beg}}, c_{\text{end}}]$ of CC start columns, for splitting one CC
  calculation into several runs which together cover all columns $0, \dots, n-1$.

#### Restricting the CC Search (Expert Options)

By default, CC grows clusters from every column $i$ using only columns $j > i$. The expert options `start` and
`cbeg`/`cend` (with `method=2` or `method=3`) and `noscan` (`method=2` only) restrict this search. Each of them
prints a `# WARNING:` to `stderr`, since the reported $d_{\min}$ is then not necessarily a lower bound on the
distance. The reported $d_{\max}$ is the weight of an actual codeword, and it is always a valid upper bound.

| Option          | CC start columns                           | Cluster growth  | Reported $d_{\min}$ valid    |
|-----------------|--------------------------------------------|-----------------|------------------------------|
| (default)       | all, $0 \le i \le n-1$                     | columns $j > i$ | always                       |
| `start=a,b,c`   | listed columns only                        | unlimited       | under a code symmetry only   |
| `cbeg`, `cend`  | $c_{\text{beg}} \le i \le c_{\text{end}}$  | columns $j > i$ | after combining split runs   |
| `noscan=1`      | all, $0 \le i \le n-1$                     | columns $j > i$ | only if `dmin=wmax` is given |

- **`start=a,b,c`** (e.g., quasi-cyclic codes): CC clusters are grown only from the listed (0-based) columns, but
  they are **not** limited to larger column indices, so that every codeword whose support contains a listed column is
  found. `start=N` is a one-element list, and `start=-1` (default) means all columns; the list is sorted, and
  duplicates are removed. The reported $d_{\min}$ (and an exact result $d_{\min} = d_{\max}$) is valid only if a code
  symmetry maps every minimum-weight codeword to a codeword containing a listed column. E.g., for a quasi-cyclic code
  with circulant blocks of size $\ell$ (columns ordered block by block), simultaneous cyclic shifts of all blocks are
  a symmetry, and it is sufficient to list one column per block, e.g., `start=0,48,96` for three blocks with
  $\ell = 48$. Cannot be combined with `cbeg`/`cend`.
- **`cbeg=[int]`, `cend=[int]`** (split runs): CC starts only from columns $c_{\text{beg}} \le i \le c_{\text{end}}$
  (default `-1`: column $0$, respectively, column $n-1$), still growing clusters using only columns $j > i$. Thus the
  $d_{\min}$ of a single run only covers codewords whose smallest column is in $[c_{\text{beg}}, c_{\text{end}}]$. To
  split a calculation, e.g., among several machines, run `cbeg=0 cend=99`, `cbeg=100 cend=199`, ..., so that every
  column $0, \dots, n-1$ is eventually listed in one of the runs, and take the minimum of $d_{\min}$ (and of the
  positive $d_{\max}$) over all runs.
- **`noscan=1`** (`method=2` only): a single CC round at $w = w_{\max}$, skipping weights $w < w_{\max}$. A codeword
  found only gives an upper bound $d_{\max}$, and $d_{\min}$ is not raised above the supplied `dmin`. The result is
  certified only if `dmin=wmax` is supplied, i.e., if all smaller weights are already known to be absent.
- **Before version 0.10.0**, `start=N` was equivalent to `cbeg=N cend=N` (clusters started from column $N$ using
  only larger columns; use `cbeg=N cend=N` to reproduce it), and `noscan=1` results were reported as certified even
  though smaller weights were not checked.
- **Python interface** (`dist_m4ri.py`): `start` takes a list of columns (`start=[0, 48]`, `start="0,48"`, or
  `start=0,48` in the CLI), while `noscan`, `cbeg`, and `cend` are disabled: they are accepted for backward
  compatibility, but ignored with a warning to `stderr` (run the `dist_m4ri` binary directly to use them). With the
  persistent cache, start-list results are stored under the separate key `<key>:start=a,b,c`; such a run is seeded
  with the bounds of the main record, and only the always valid $d_{\max}$ and codewords are merged into the main
  record. The expert option `trust_start` (CLI: `--trust-start` or `trust_start=1`) accepts the start-list results as
  valid: $d_{\min}$ and `rw_steps` are merged as well, and an existing exact start-list record is copied to the main
  record (and the cache file is saved) even if no calculation is needed. For DEM input, the start columns index the
  error mechanisms in their order in the flattened DEM (after the `pmin` cutoff). For CSS codes, the same list is
  used in both sectors.
- **Possible optimization** (not implemented): clusters grown from a listed column could skip the columns listed
  before it, since all codewords containing those columns have already been enumerated.

### 3. Bracketing Mode (`method=3`, default)
Dynamically partitions the available thread pool between CC (pushing $d_{\min}$ up) and RW (pulling $d_{\max}$ down) to
determine the exact code distance as quickly as possible.

#### Dynamic Thread Allocation & Role of `dexp`:
1. **Target Search Depth**:
   - Before RW discovers a candidate codeword, $d_{\max}$ is unknown. Providing `dexp=D` (alias: `dest=D`) tells the
     coordinator to plan CC verification up to target weight $D - 1$ (or $D$).
2. **Predictive Workload Modeling**:
   - The coordinator measures empirical single-thread speed for RW steps ($t_{\text{RW}}$) and fits exponential growth
     to completed CC rounds to estimate time $T_{\text{CC}}$ required to reach $\min(d_{\max}, d_{\exp})$.
   - It computes the remaining work ratio:
     $$\text{ratio} = \frac{T_{\text{CC}}}{T_{\text{CC}} + T_{\text{RW}}}$$
     and dynamically splits threads at each round:
     $$N_{\text{CC}} = \text{round}(N_{\text{threads}} \times \text{ratio}),$$
     $$N_{\text{RW}} = N_{\text{threads}} - N_{\text{CC}}$$
     where $N_{\text{RW}}$ cannot exceed the number of threads allowed to run RW (see
     [Multithreading, Throttling & Batch Sizing](#4-multithreading-throttling--batch-sizing)).
   - **Small $d_{\exp}$**: CC requires little work, so only 1–2 threads run CC while the majority maximize RW sampling
     speed.
   - **Heavier rounds**: As $w$ grows toward $d_{\exp}$, $T_{\text{CC}}$ increases and additional threads are shifted
     to CC to ensure both algorithms converge on the exact distance simultaneously.
3. **Adaptive Early Cutoff & RW Convergence**:
   - When RW satisfies the `min_hits` convergence criterion, RW workers stop early and yield 100% of threads to CC
     to finish certifying $d_{\min}$.
   - Likewise, once all RW `steps` have been claimed, RW threads join the current CC round instead of waiting for it
     to finish. The CC work of each round is measured as the total CC thread time.
   - If CC reaches $w > d_{\exp}$ before a codeword is found, CC halts and yields 100% of threads to RW.
   - If a projected CC round is estimated to exceed the remaining `timeout`, the coordinator terminates CC early and
     devotes remaining time entirely to RW.

Relevant parameters:
- `dexp=[int]` (alias: `dest=[int]`): Expected code distance to guide target search depth and thread allocation.
- `threads=[int]`: Maximum number of worker threads (default: hardware concurrency; subject to throttling unless
  `nothrottle=1` is specified).
- `nothrottle=[int]`: Disable automatic thread throttling (default: 0; CLI flag: `--no-throttle`).
- `chunk_size=[int]` (alias: `batch=[int]`): RW step batch chunk size (default: 0 for automatic adaptive sizing).
  When `chunk_size=0` (default in Python and binary), chunk size is chosen adaptively based on $n$ and step budget:
  - Small codes ($n < 500$): 50 steps/chunk (250 for $\ge 10^3$ steps; 500 for $\ge 5 \times 10^4$ steps).
  - Medium codes ($500 \le n < 5000$): 50 steps/chunk (100 for $\ge 10^4$ steps).
  - Large matrices ($n \ge 5000$): 25 steps/chunk (50 for $\ge 10^4$ steps;
    bounds timeout overshoot to $\le 2\text{--}3$s).
- `timeout=[sec]`: Maximum execution time in seconds (default: 60.0; set to `0` for infinite / no timeout).
- `steps=[int]`: Maximum total RW steps (default: 100000; set to `0` to run pure CC via bracketing coordinator).
- `min_hits=[int]`: Minimum hit count per minimum-weight codeword for early RW termination (default: 5; 0 to disable).
- `dW=[int]`: Extra weight window above $d_{\min}$ to continue collecting codewords ($w \le d_{\min} + \text{dW}$).

### 4. Multithreading, Throttling & Batch Sizing
The parameter `threads=[int]` specifies the **maximum** number of worker threads to allocate.
By default (`nothrottle=0`), automatic heuristics clamp thread usage to avoid overhead and resource thrashing:
- **Small-Code Throttling**: When $n < 100$ or $r \cdot n < 100,000$, threads are automatically clamped to $\le 4$
  (and $\le 16$ for $n < 300$), eliminating thread creation and lock contention overhead.
- **Large-Matrix Memory Throttling (RW threads only)**: For massive matrices (e.g. circuit DEMs with
  $n > 20,000, r > 5,000$), the number of RW threads is automatically limited so that their total dense working memory
  (copies of $H$ and $H^T$ per RW thread, or the `ksub` $\times n$ subspace matrix with `ksub>0`) stays under ~1.5 GB,
  avoiding DRAM bus and CPU cache thrashing. CC needs no dense matrices: `method=2` is not limited, and in `method=3`
  only the RW share is limited, while CC rounds can use all threads.
- **Small Step Counts (RW threads only)**: At most $\lceil \text{steps} / 10 \rceil$ threads run RW (all threads in
  `method=1`, the RW share in `method=3`, where CC rounds can still use all threads; `steps=0` in `method=3` runs pure
  CC without RW threads).
- **Adaptive RW Chunk Sizing**: Setting `chunk_size=0` (or omitting it) enables adaptive batching, reducing atomic CAS
  contention by up to 10x during extended searches ($10^4$ or $10^5$ steps).
- **Thread Starvation Prevention**: If `chunk_size` exceeds $\lceil \text{steps} / N_{\text{RW}} \rceil$, where
  $N_{\text{RW}}$ is the number of RW threads, the chunk is automatically clamped so that a single thread cannot
  monopolize all steps, ensuring all cores run concurrently.
- **CSS Codeword Suffixing**: In CSS mode, specifying `outC="cws.nz"` automatically saves $X$-codewords to `cws_X.nz`
  and $Z$-codewords to `cws_Z.nz` (preventing mixed sectors in a single file). Specifying `finC="cws.nz"` automatically
  resolves `cws_X.nz` and `cws_Z.nz` (or separates mixed files in-flight).
- **Manual Override**: Pass `nothrottle=1` (or `--no-throttle` in Python) to force `dist_m4ri` to use the exact number
  of requested threads without throttling, and specify `chunk_size=N` to override batch sizing.

---

## Confinement Profile

By default, `smax=0` (confinement tracking disabled for maximum CC speed; in `.mtx` mode with `debug&1`, a notice is
logged to `stderr`). When `smax > 0` (e.g. `smax=5`), the CC algorithm tracks the minimum non-zero syndrome weight
observed for each cluster weight $w$:

```bash
$ ./src/dist_m4ri method=2 finH=./examples/surf_d5_H.mmx finL=./examples/surf_d5_L.mmx wmax=4 smax=5 debug=0 threads=4
# confinement: 1,1,1,1
5 0 0
```

With `debug=1`, detailed per-weight lines are printed to `stderr`:
```text
# w=1 min non-zero syndrome weight 1
# w=2 min non-zero syndrome weight 1
# w=3 min non-zero syndrome weight 1
# w=4 min non-zero syndrome weight 1
```

---

## Codeword Export (`outC` / `finC`)

- **`outC=[file.nz]`**: Saves all unique discovered codewords in standard **NZLIST** format:
  ```text
  %% NZLIST
  % generated by dist_m4ri
  <weight> <col_1> <col_2> ... <col_weight>
  ```
  *(Indices are 1-based).*
- **`finC=[file.nz]`**: Reads initial codewords from a file to initialize $d_{\max}$ and the codeword hash table.
- **`dW=[int]`**: When set (e.g. `dW=1`), preserves and exports codewords of weight up to $w \le d_{\min} + \text{dW}$.
- **`maxC=[int]`**: Limits collection to at most `maxC` unique codewords.

---

## Command-Line Usage and Help System

`dist_m4ri` provides a three-tier help system:

1. **Short help on error**: Printed to `stderr` when invoked without arguments or with an unrecognized parameter,
   listing all allowed parameters concisely.
2. **Standard `--help` (fits 80 rows)**: Displays the most commonly used parameters with clear descriptions, bunches
   extra available parameters at the bottom, and points to `--morehelp`.
3. **Full `--morehelp`**: Displays complete, detailed descriptions of all command-line parameters, debug bitmask flags,
   and output formats.

### Standard Help (`dist_m4ri --help`)

```text
$ ./src/dist_m4ri --help
./src/dist_m4ri (version 0.10.1): calculate distance of a classical or quantum CSS code
Usage: ./src/dist_m4ri [method=1|2|3] [parameter=value ...]

Calculation method:
  method=[int]       1: Random Window (RW) algorithm (upper bound)
                     2: Connected Cluster (CC) algorithm (lower bound / exact)
                     3: Bracketing mode (concurrent RW and CC) (default: 3)

Input matrices (Matrix Market .mmx/.mtx format or Stim DEM):
  finH=[file]        Parity check matrix H (classical) or Hx (CSS quantum)
  finG=[file]        Hz check matrix (quantum CSS code only)
  finL=[file]        Lx logical operator matrix (quantum CSS code only)
                     Note: Either L=Lx or G=Hz is required for quantum CSS codes
  fin=[str]          Base name for CSS matrices (loads ${fin}X.mtx, ${fin}Z.mtx)
  fdem=[file]        Stim detector error model (DEM) file
  pmin=[float]       Minimum error probability threshold to keep for DEM (0.0)
  classical=[0|1]    1: classical code (Hx only), 0: quantum CSS (auto-detected)

Distance bounds and guidance:
  dmin=[int]         Certified lower bound on distance (CC starts from dmin) (1)
  dmax=[int]         Known upper bound on distance (RW ignores cw wt >= dmax) (0)
  dexp=[int]         Expected distance for method=3 thread allocation (alias: dest) (0)

Search limits and stopping criteria:
  steps=[int]        Maximum RW decoding steps / information sets (100000)
  wmax=[int]         Maximum cluster weight to search in CC (0=until bound/timeout)
  wmin=[int]         Stop immediately if cw with weight <= wmin is found (1)
  min_hits=[int]     Stop RW when min-wt cws (at least cov_cws) hit >= min_hits (5)
  timeout=[sec]      Execution timeout in seconds, 0 for infinite (60.0)

Multithreading and RW optimization:
  threads=[int]      Max worker threads to use (0: auto CPU count) (0)
  ksub=[int]         Subspace dimension sampled from ker(H) for RW (0: full H) (0)
  kwin=[int]         Localized column permutation window size W (0: auto/hybrid) (0)

Codeword collection:
  outC=[file]        Export found minimum-weight codewords to file (.nz format)
  finC=[file]        Import initial codewords from file (.nz format)
  maxC=[int]         Maximum number of codewords to collect (0 for unlimited) (0)
  dW=[int]           Collect codewords up to weight dmin + dW (default: 0)

Extra parameters (see --morehelp for details):
  smax=[int] (0)         Max syndrome weight for confinement profile (0 to disable)
  noscan=[0|1] (0)       Expert, method 2: CC at w=wmax only (no lower bound!)
  start=[list] (-1)      Expert: CC from listed columns only, e.g. start=0,48
  cbeg/cend=[int] (-1)   Expert: CC start column range, for split runs
  nothrottle=[0|1] (0)   Disable thread throttling (also --no-throttle)
  chunk_size=[int] (0)   RW batch chunk size (0: auto, alias: batch)
  win_mode=[0|1] (0)     Window mode: 0=Tanner BFS, 1=index proximity
  cov_cws=[int] (100)    Max min-wt cws tracked in hash for min_hits stop
  refresh=[int] (0)      RW steps between adaptive ker(H) basis refreshes (auto 5000 if ksub>0)
  seed=[int] (0)         RNG seed [0 for time(NULL)]
  debug=[int] (3)        Debug bitmask (0: silent, 1: general, 2: verbose, ...)

Help options:
  -h, --help         Display this help message (commonly used parameters)
  --morehelp         Display full help with all parameter descriptions
  --version          Display program version
```

### Full Parameter Listing (`dist_m4ri --morehelp`)

Use `./src/dist_m4ri --morehelp` to view exhaustive parameter explanations, including all `debug` bitmask flags
(e.g., `debug=0` for silent mode, `debug=1` for general info, `debug=2` for thread timing, `debug=16` for new
codewords, `debug=32` for matrix dumps).

### CLI Examples

```bash
# 1. Classical linear code using 8 threads in bracketing mode (method=3 is default)
$ ./src/dist_m4ri finH=./examples/c204H.mmx dest=10 steps=100000 threads=8 debug=0
8 8 150

# 2. Stim Detector Error Model (DEM) with timeout and codeword export
$ ./src/dist_m4ri fdem=./examples/surf_d3.dem dexp=3 outC=cws.nz threads=4 debug=0
3 3 50

# 3. Quantum CSS code (Hx and Lx) using pure CC search up to wmax=5 with confinement
$ ./src/dist_m4ri method=2 finH=./examples/surf_d5_H.mmx finL=./examples/surf_d5_L.mmx wmax=5 smax=5 debug=0 threads=4
# confinement: 1,1,1,1,1
5 5 0
```

---

## Python Wrapper (`dist_m4ri.py`)

A Python module [`dist_m4ri.py`](dist_m4ri.py) is included for high-level scripting, NumPy/SciPy integration, and Stim
interoperability without manual threading overhead.

### Key Python Functions

- `compute_classical_distance(H, ...)`: Minimum distance of a classical linear code (from NumPy 2D array, SciPy sparse
  matrix, or `.mtx` file).
- `compute_css_distance(Hx, Hz, Lx=None, Lz=None, ...)`: Distance $d = \min(d_X, d_Z)$ of a CSS quantum code.
- `compute_dem_distance(dem=None, circuit=None, simple=None, full=False, basis=None, rounds=None, out_dir=None,`:
  `out_dem=None, out_stim=None, ...)`:
  Minimum distance directly from a `stim.DetectorErrorModel`, `stim.Circuit`, `.dem` file, or `.stim` circuit file:
  - Always performs **ancilla and data basis tracking** (`classify_qubits_thorough`) across virtual SWAP permutations
    and Clifford gates to identify data qubits, $X$-sector and $Z$-sector ancillas, routing qubits, and local data basis
    rotations (automatically recognizing both standard CSS codes and locally basis-rotated CSS codes such as **XZZX**).
  - **`--simple` vs `--full`**: For CSS and locally rotated CSS `.stim` circuits, `--simple` (`simple=True`) is enabled
    by default, stripping minority-basis detectors (`strip_minority_detectors`) to keep only primary-basis detectors and
    produce a smaller, simpler DEM. Pass `--full` (`full=True`) to retain all detectors.
  - **`--rounds [N]`**: Constructs the DEM with `N` repetitions in the circuit's `REPEAT` block (`rounds=2` by default
    when `--rounds` is passed without a number). If the `.stim` circuit does not contain a `REPEAT` block, issues a
    warning to `stderr` and continues with the original circuit.
  - **`--out-dir DIR`, `--out-dem [FILE]`, `--out-stim [FILE]`**: Optionally saves the constructed DEM and/or the
    processed (noisy, round-adjusted, detector-filtered) `.stim` circuit. When `--out-dem` or `--out-stim` is passed
    without a filename (or `out_dem=True` / `out_stim=True`), filenames default to `<basename>_simp.dem` /
    `<basename>_full.dem` and `<basename>_simp.stim` / `<basename>_full.stim`. By default, `--out-dir` defaults to the
    directory of the input file.
- `classify_qubits_thorough(circuit, ...)` / `strip_minority_detectors(circuit, basis, ...)` /
  `set_circuit_rounds(circuit, rounds)`: Helpers for Stim circuit Pauli basis tracking, minority-detector stripping, and
  `REPEAT` block round adjustment.
- `has_noise(circuit)` / `add_noise(circuit, noise_prob=0.001)`: Inspects a `stim.Circuit` for noise instructions and
  injects phenomenological `DEPOLARIZE1` / `DEPOLARIZE2` / reset-flip / measurement-flip noise into noiseless circuits.
- `read_sparse_vectors(filepath)`: Parses NZLIST files into lists of 0-based integer support indices.
- Distance caching: `enable_distance_cache()`, `disable_distance_cache()`, `clear_distance_cache()`,
  `get_cached_distance(..., start=None)`.
- Expert CC options of all `compute_*_distance()` functions: `start` (list of CC start columns, e.g. `start=[0, 48]`;
  separate cache record `<key>:start=a,b,c`) and `trust_start=True` (CLI: `--trust-start`; accept the start-list
  results as valid and copy them to the main cache record). The binary options `noscan`, `cbeg`, and `cend` are
  disabled in Python (ignored with a warning to `stderr`); see
  [Restricting the CC Search](#restricting-the-cc-search-expert-options).
- Optional solver backend: `solver="codedistance"` (uses the `codedistance` library if installed).

### Python Example

```python
import numpy as np
import stim
import dist_m4ri

# 1. Classical Code Distance
H = np.array([
    [1, 0, 0, 1, 1, 0, 1],
    [0, 1, 0, 1, 0, 1, 1],
    [0, 0, 1, 0, 1, 1, 1]
], dtype=np.int8)
d = dist_m4ri.compute_classical_distance(H, threads=4)
print(f"Hamming code distance: {d}")  # 3

# 2. Stim DEM / Circuit Distance with Codewords
circuit = stim.Circuit.generated(
    "surface_code:rotated_memory_z",
    rounds=3,
    distance=3,
    after_clifford_depolarization=0.001
)
dist, dist_list, cws = dist_m4ri.compute_dem_distance(circuit=circuit, do_cws=True, threads=4)
print(f"Surface code distance: {dist}, found {len(cws)} minimum-weight error mechanisms")

# 3. CSS Quantum Code
dist, d_x, d_z = dist_m4ri.compute_css_distance(
    Hx="examples/surf_d5_H.mmx",
    Hz="examples/surf_d5_H.mmx",
    Lz="examples/surf_d5_L.mmx",
    Lx="examples/surf_d5_L.mmx",
    d_exp=5,
    threads=4
)
print(f"CSS distance: {dist}")  # 5
```

---

## Compilation & Testing

### Prerequisites
- Recent `gcc` with POSIX threads support (`-pthread`).
- `libm4ri-dev` linear algebra library:
  ```bash
  sudo apt-get update -y
  sudo apt-get install -y libm4ri-dev
  ```

### Build Targets

```bash
cd src

# Compile both multithreaded dist_m4ri and single-threaded dist_m4ri_old
make all

# Run full C test suite (74 tests)
make test
```

### Python Unit Tests

```bash
pytest tests/test_dist_m4ri.py -v
```

### Benchmark Suite (`benchmark/`)

The [`benchmark/`](benchmark/BENCHMARK.md) directory provides quantum CSS codes ($d \in [6, 36]$) and Stim circuits
(Gross code, Mitten codes, honeycomb/color codes, bivariate bicycle codes) along with instructions for generating full
and stripped DEMs via `dist_m4ri.add_noise`.

---

## References

If you use this program, please cite:

*   A. Dumer, A. A. Kovalev, and L. P. Pryadko, "Distance verification for classical and quantum LDPC codes,"
    *IEEE Transactions on Information Theory*, vol. 63, no. 7, pp. 4675-4690, 2017.
    [doi:10.1109/TIT.2017.2690381](https://doi.org/10.1109/TIT.2017.2690381).

Other related papers and software:

*   **vecdec Repository** (Random Information Set (RW) algorithm with error weights/probabilities):
    [QEC-pages/vecdec](https://github.com/QEC-pages/vecdec).

*   **QDistRnd GAP Package** (Random Information Set algorithm for quantum codes over arbitrary finite fields):
    L. P. Pryadko, V. A. Shabashov, and V. K. Kozin, "QDistRnd: A GAP package for computing the distance of quantum
    error-correcting codes," *Journal of Open Source Software*, vol. 7, no. 71, p. 4120, 2022.
    [doi:10.21105/joss.04120](https://doi.org/10.21105/joss.04120).

*   **Performance Comparison**:
    M. Webster, A. Jacob, and O. Higgott, "Distance-Finding Algorithms for Quantum Codes and Circuits,"
    arXiv:2603.22532 [quant-ph], 2026. [arXiv:2603.22532](https://arxiv.org/abs/2603.22532).


