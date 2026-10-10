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
  estimate (`dexp`/`dest`), remaining RW steps and their measured time, timeout, and the measured scaling
  characteristics of CC.

For a classical binary linear code, only the parity-check matrix $H$ is needed.

For a quantum CSS code, matrix $H_X$ and either $H_Z$ (as `finG`) or logical operators $L_X$ (as `finL`) are needed.

Alternatively, a detector error model file from `stim` can be specified using `fdem=[str]`.

All matrices with entries in $\text{GF}(2)$ have $n$ columns and obey the orthogonality conditions:
$$H_X H_Z^T = 0,\quad H_X L_Z^T = 0,\quad L_X H_Z^T = 0,\quad L_X L_Z^T = I.$$
The last condition is not required: it is sufficient that $L_X$, $L_Z$, and $L_X L_Z^T$ have the same full row rank
$k = n - \mathrm{rank}\,H_X - \mathrm{rank}\,H_Z$, the number of encoded qubits (printed with `debug=16`).

---

## Output Format & Stream Separation

`dist-m4ri` strictly separates machine-parseable results from informational progress logs:

- **`stdout`**: Outputs three space-separated integers:
  ```text
  dmin dmax rw_steps
  ```
  - `dmin - 1` is the maximum cluster size analyzed without success by CC (`dmin = dmax` if CC found a minimum-weight
    codeword).
  - `dmax` is the weight of the smallest non-trivial codeword found (by RW or CC, or read with `finC`), or the supplied
    upper bound `dmax=[int]` if smaller (`0` if none).
  - `rw_steps` is the number of completed RW steps across all threads (`0` if CC found a minimum-weight codeword, or
    if RW did not run in `method=2`).
  - When `dmin = dmax = d`, the exact code distance is confirmed.
  - With a stop target `dstop=U` (`method=2` or `3`), the run ends once CC has certified `dmin >= U`, unless a codeword
    of weight `< U` is found; `dmax` is then the weight of the lightest codeword found (`0` if none), never `U`, since
    `dstop` is not an upper bound (see [Bracketing Mode](#3-bracketing-mode-method3-default)).

> **Note on Compatibility**: This 3-number output format (`dmin dmax rw_steps`) is specific to the multithreaded
> `dist_m4ri` and is incompatible with the legacy single-threaded `dist_m4ri_old` (which returned a single integer `d`
> or `-w`).

- **`stderr`**: Receives all status banners, thread balancing reports, CC round timings, RW discovery logs, warnings,
  and confinement profiles, as selected by the `debug` bitmap (see [Debug Output](#debug-output-debugint)).

---

## How the Methods Work

### 1. Multithreaded RW Algorithm (`method=1`)
Searches for low-weight non-trivial binary codewords $c$ such that $Hc = 0$ and $Lc \neq 0$.
Threads independently generate random column permutations, compute Gaussian elimination to find information sets,
and extract candidate dual-row codewords. Only the candidates lighter than the current weight limit are formed (the
weight of each candidate is counted first), and they are checked for $Lc \neq 0$ with bit masks of the columns of $L$,
in $O(|c|)$ operations: for DEMs, most such candidates are trivial. When a lighter codeword is discovered, all worker
threads atomically update the global upper bound $d_{\max}$. The found codewords are kept in a hash table with their
hit counts (the number of times RW found each codeword): the collection window for `outC` (weight up to
$w_{\min} + \text{dW}$, where $w_{\min}$ is the minimum weight found) and, for `min_hits`, a representative set of
`cov_cws` lowest-weight codewords, where heavier codewords are gradually replaced as lighter ones are found.

Relevant parameters:
- `steps=[int]`: Total number of information sets / RW rounds across all threads (default: 100000).
- `min_hits=[int]`: QDistRnd-style automatic convergence stopping criterion (default: 5; set 0 to disable). Stops RW
  when the average number of hits per codeword $\langle n\rangle$ reaches `min_hits`, both for the codewords of the
  minimum weight found and for the representative set of `cov_cws` lowest-weight codewords (heavier codewords make the
  estimate more conservative). A single codeword suffices (e.g., a classical code with $k=1$). A lighter codeword is
  then missed with probability about $e^{-\langle n\rangle}$; see
  [RW Convergence](#5-rw-convergence-empirical-test-and-probability-to-find-a-codeword) for this empirical estimate,
  its test (a warning at the end with `debug&1`), and the information-set estimate.
- `cov_cws=[int]`: Size of the representative set of lowest-weight codewords for `min_hits` (default: 100; 0 for the
  minimum-weight codewords only). Without `outC` and `maxC`, at most `cov_cws` minimum-weight codewords are tracked
  (all with `cov_cws=0`). For codes with fewer than `cov_cws` minimum-weight codewords, a smaller `cov_cws` may stop
  RW earlier.
- `kwin=[int]` (alias: `win=[int]`): Localized column permutation window size $W$: the random column order starts with
  $W$ columns near a random column, which are thus preferred as pivots (default: 0; for $n \ge 500$, every other RW
  step then uses $W = \min(512, 3n/4)$ consecutive columns, as with `win_mode=1`, and the other steps use uniform
  random permutations).
- `win_mode=[0|1]`: Window construction mode when `kwin > 0`: `0` for Tanner graph BFS neighborhood, `1` for contiguous
  column index window (default: 0).
- `ksub=[int]`: **Experimental, should not be used** (a warning is printed whenever `ksub>0` is given for RW, with or
  without `min_hits`). Subspace dimension sampled from $\ker(H)$ for cache-resident RW elimination (default: 0 for
  full-matrix RW; automatically falls back to `ksub=0` with a warning when $m < \nu = \dim\ker(H)$, e.g., for DEMs
  and most quantum codes). Each step reduces the span of `ksub` distinct rows of a fixed basis of $\ker(H)$, so that a
  few codewords are found much more often than others: then $e^{-\langle n\rangle}$ is optimistic, and `min_hits` may
  end RW early.
- `refresh=[int]`: RW step interval for adaptive $\ker(H)$ basis refresh via low-weight codeword exchange and
  re-echelonization, only with the experimental `ksub>0` (default: 5000 when `ksub > 0`, 0 to disable).
- `wmin=[int]`: Minimum distance of interest (stop immediately when a codeword of weight $w \le w_{\min}$ is found;
  not when collecting codewords with `outC` or `maxC`).
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
soon as $d_{\min} = d_{\max}$, unless codewords are collected (`outC` or `maxC`). A stop target `dstop=U` ends the scan
once $d_{\min} \ge U$ (with `outC`, after the rounds $w = U, \dots, U + \text{dW}$), and $d_{\max}$ remains the weight
of the lightest codeword found (`0` if none).
With a `timeout`, each round is started even if it is predicted not to finish in time (its CC work, the total CC
thread time, is extrapolated from the last two rounds with their growth factor, clamped to $[2, 10]$; a note is
printed): a round which cannot be completed may still find a codeword of weight $w = d_{\min}$, i.e., the exact
distance.

Relevant parameters:
- `wmax=[int]`: Maximum cluster weight to search (optional with `timeout>0`, the default, with `dmax>0`, where the
  rounds end once $d_{\min} = d_{\max}$, or with `dstop>0`, where they end once $d_{\min} \ge$ `dstop`; otherwise
  required for `method=2`).
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
   - Once a codeword has been found, CC rounds continue up to $w = d_{\max} - 1$: the round $w = d_{\max} - 1$
     certifies $d_{\min} = d_{\max}$ (with `outC`, CC continues with the rounds $w = d_{\max}, \dots, d_{\max} +
     \text{dW}$, which enumerate all codewords to export). Before that, CC is limited only by `wmax` and the `timeout`.
   - `dexp=D` (alias: `dest=D`) is a hint: as long as RW has not found any codeword, CC rounds $w > D$ run only on the
     threads which cannot run RW (see
     [Multithreading, Throttling & Batch Sizing](#4-multithreading-throttling--batch-sizing)), i.e., CC pauses if RW
     can use all threads. CC resumes as soon as RW finds a codeword, or when RW ends.
   - A stop target `dstop=U` (not an upper bound) takes the place of $d_{\max}$ as long as no codeword of weight
     $< U$ is known: the run ends once the round $w = U - 1$ certifies $d_{\min} \ge U$ (with `outC`, after the rounds
     $w = U, \dots, U + \text{dW}$), without waiting for RW, and $d_{\max}$ is the weight of the lightest codeword found
     (`0` if none). E.g., `dstop=U wmin=U-1` decides whether $d \ge U$, and for a CSS code with
     $d = \min(d_X, d_Z) \le U$ known (e.g., from the other sector), `dstop=U` is all that is needed.
2. **Predictive Workload Modeling**:
   - The RW threads measure the RW step time $t_{\text{RW}}$ (thread time per step) continuously. The CC work of each
     round is measured as the total CC thread time, and the work of the next round is extrapolated with the growth
     factor of the last two rounds (clamped to $[2, 10]$).
   - For each round, the coordinator compares the remaining work in thread-seconds: $T_{\text{CC}}$ of the rounds
     $w, w+1, w+2$ (up to $d_{\max} - 1$, or up to $d_{\exp}$ before a codeword is found), and
     $T_{\text{RW}} = t_{\text{RW}} \times$ (remaining RW steps, at most 2000 once $d_{\max}$ is known), but no more
     than the RW threads can do before the `timeout`. It splits the threads as
     $$\text{ratio} = \frac{T_{\text{CC}}}{T_{\text{CC}} + T_{\text{RW}}},$$
     $$N_{\text{CC}} = \text{round}(N_{\text{threads}} \times \text{ratio}),$$
     $$N_{\text{RW}} = N_{\text{threads}} - N_{\text{CC}}$$
     with $1 \le N_{\text{CC}} \le N_{\text{threads}} - 1$, where $N_{\text{RW}}$ cannot exceed the number of threads
     allowed to run RW (see [Multithreading, Throttling & Batch Sizing](#4-multithreading-throttling--batch-sizing)).
   - **Small rounds** (estimated CC work below 5 ms): only 1–2 threads run CC while the majority maximize RW sampling
     speed.
   - **Heavier rounds**: as the CC work grows, additional threads are shifted to CC, so that both algorithms converge
     on the exact distance simultaneously.
   - **The certifying round** $w = d_{\max} - 1$ (and any later codeword collection round) runs on all threads, since
     RW can no longer improve the result.
   - **The first round at a supplied `dmin`** (no measured round yet) runs on half of the threads.
   - **During a round**, CC threads are added when RW finds a new $d_{\max}$ (e.g., if the round now certifies
     $d_{\min} = d_{\max}$), and all threads join a round which takes more than 1.5 times its predicted time. A busy
     RW thread joins a CC round when its current RW batch is done.
3. **Timeout, Early Cutoff & RW Convergence**:
   - With a `timeout`, a CC round gets enough CC threads to finish within 3/4 of the remaining time. A round predicted
     not to finish in time on all threads is not started while RW runs; once RW has ended, the run ends instead. With
     `steps=0` (no RW), as in `method=2`, such a round is started anyway, since it may still find a codeword of weight
     $w = d_{\min}$.
   - The coordinator never ends CC while RW runs. Whenever no CC round can be started (CC done up to `wmax`, a round
     predicted to exceed the `timeout`, or CC paused by `dexp`), it waits for RW and re-plans when RW finds a new
     $d_{\max}$ or ends. The run ends once RW has ended and no CC round can be started.
   - When RW satisfies the `min_hits` convergence criterion, RW workers stop early and yield 100% of threads to CC
     to finish certifying $d_{\min}$.
   - Likewise, once all RW `steps` have been claimed, RW threads join the current CC round instead of waiting for it
     to finish.

Relevant parameters:
- `dexp=[int]` (alias: `dest=[int]`): Expected code distance, a hint for the thread allocation before RW finds a
  codeword (see above).
- `dstop=[int]`: Stop target for the lower bound (default: 0, off; see above, also for `method=2`). Unlike `dmax`, it is
  not an upper bound, and it is never reported as `dmax`. A supplied `dmin >= dstop` ends the run without a search
  (unless codewords are collected with `outC` or `maxC`); otherwise, `dstop` has no effect with `method=1`.
- `threads=[int]`: Maximum number of worker threads (default: hardware concurrency; subject to throttling unless
  `nothrottle=1` is specified).
- `nothrottle=[int]`: Disable automatic thread throttling (default: 0; CLI flag: `--no-throttle`).
- `chunk_size=[int]` (alias: `batch=[int]`): RW step batch chunk size (default: 0 for automatic adaptive sizing).
  When `chunk_size=0` (default in Python and binary), chunk size is chosen adaptively based on $n$ and step budget:
  - Small codes ($n < 500$): 50 steps/chunk (250 for $\ge 10^3$ steps; 500 for $\ge 5 \times 10^4$ steps).
  - Medium codes ($500 \le n < 5000$): 50 steps/chunk (100 for $\ge 10^4$ steps).
  - Large matrices ($n \ge 5000$): 25 steps/chunk (50 for $\ge 10^4$ steps).
- `timeout=[sec]`: Maximum execution time in seconds (default: 60.0; set to `0` for infinite / no timeout). RW threads
  check the timeout (and the other stop conditions) between RW steps, and also within a step (every 16 columns of the
  Gaussian elimination), so that even on huge matrices, where a single RW step may take many seconds, the run ends
  shortly after the timeout.
- `steps=[int]`: Maximum total RW steps (default: 100000; set to `0` to run pure CC via bracketing coordinator).
- `min_hits=[int]`: Average number of hits per lowest-weight codeword for early RW termination (default: 5; 0 to
  disable; see [method=1](#1-multithreaded-rw-algorithm-method1)).
- `dW=[int]`: Extra weight window above $d_{\min}$ to continue collecting codewords ($w \le d_{\min} + \text{dW}$).

### 4. Multithreading, Throttling & Batch Sizing
The parameter `threads=[int]` specifies the **maximum** number of worker threads to allocate.
By default (`nothrottle=0`), automatic heuristics clamp thread usage to avoid overhead and resource thrashing:
- **Small-Code Throttling** (all threads), to avoid the thread creation and lock contention overhead on tiny
  workloads; with $m$ the number of rows of $H$ and $W_{\text{RW}} = m \cdot n \cdot \text{steps}$ the RW workload
  ($W_{\text{RW}} = 0$ in `method=2`): at most 4 threads for $n < 60$, or for $m \cdot n < 2 \cdot 10^4$ with
  $W_{\text{RW}} < 5 \cdot 10^7$; otherwise at most 16 threads for $n < 150$ or $m \cdot n < 10^5$, if
  $W_{\text{RW}} < 2 \cdot 10^8$.
- **Large-Matrix Memory Throttling (RW threads only)**: For massive matrices (e.g. circuit DEMs with
  $n > 20,000, m > 5,000$), the number of RW threads is automatically limited so that their total dense working memory
  (copies of $H$ and $H^T$ per RW thread, or the `ksub` $\times n$ subspace matrix with `ksub>0`) stays under ~1.5 GB,
  avoiding DRAM bus and CPU cache thrashing: above 15 MB per RW thread (about $m \cdot n > 6 \cdot 10^7$ for
  full-matrix RW), at most 1.5 GB divided by the memory per thread, but at least 2 and at most 32 threads run RW. CC
  needs no dense matrices: `method=2` is not limited, and in `method=3` only the RW share is limited, while CC rounds
  can use all threads.
- **Small Step Counts (RW threads only)**: At most $\lceil \text{steps} / 10 \rceil$ threads run RW (all threads in
  `method=1`, the RW share in `method=3`, where CC rounds can still use all threads; `steps=0` in `method=3` runs pure
  CC without RW threads).
- **Adaptive RW Chunk Sizing**: Setting `chunk_size=0` (or omitting it) enables adaptive batching, reducing atomic CAS
  contention by up to 10x during extended searches ($10^4$ or $10^5$ steps).
- **Thread Starvation Prevention**: If `chunk_size` exceeds $\lceil \text{steps} / N_{\text{RW}} \rceil$, where
  $N_{\text{RW}}$ is the number of RW threads, the chunk is automatically clamped so that a single thread cannot
  monopolize all steps, ensuring all cores run concurrently.
- **Manual Override**: Pass `nothrottle=1` (or `--no-throttle` in Python) to force `dist_m4ri` to use the exact number
  of requested threads without throttling, and specify `chunk_size=N` to override batch sizing.

With `debug&8`, the binary reports the throttling and the thread allocation of every CC round, and with `debug&4` a
periodic status line shows the progress of the current CC round. A CC round in `method=3` may have to wait until the RW
threads assigned to it finish their current RW batch.

### 5. RW Convergence: Empirical Test and Probability to Find a Codeword

The RW upper bound $d_{\max}$ is not certified unless CC confirms it. Two estimates of the probability
$P_{\text{fail}}$ that RW missed a lighter codeword are available.

**Empirical estimate and test**, as in the
[QDistRnd manual, Sec. 3.3](https://qec-pages.github.io/QDistRnd/doc/chap3.html#X7DA7BD6F7A61E553). The probability
that a RW step finds a given codeword is expected to depend only on its weight, and to decrease with the weight. If,
after $N$ steps, the $m$ distinct codewords of the minimum weight $w$ have been found $n_1, \dots, n_m$ times, the
probability per step is estimated as $\lambda_w = \sum_i n_i / (N m)$, and a codeword of a smaller weight is missed
with probability

$$
P_{\text{fail}} < (1 - \lambda_w)^N < e^{-N\lambda_w} = e^{-\langle n\rangle}, \qquad
\langle n\rangle = \frac{1}{m} \sum_{i=1}^m n_i .
$$

RW stops once $\langle n\rangle$ reaches `min_hits`: $e^{-5} \approx 0.7\%$ for the default `min_hits=5`, and
$e^{-10} \approx 4.5 \cdot 10^{-5}$ for `min_hits=10`. Here $\langle n\rangle$ is the smaller of the averages over the
minimum-weight codewords and over the representative set of `cov_cws` lowest-weight codewords (heavier codewords only
make the estimate more conservative). The assumption of equal probabilities for the codewords of the same weight can
be tested: QDistRnd uses Pearson's statistic $X^2 = (m / n_{\text{tot}}) \sum_i n_i^2 - n_{\text{tot}}$, where
$n_{\text{tot}} = \sum_i n_i$, which is distributed as $\chi^2_{m-1}$ for large $\langle n\rangle$. This is the
index-of-dispersion test of the hit counts against the Poisson distribution. `dist_m4ri` applies it to each weight
class of the representative set, using the zero-truncated Poisson distribution (codewords never found are not known),
and with `debug&1` warns at the end if the hit counts are strongly non-uniform: a p-value below $10^{-3}$, a relative
spread of the hit rates of at least 1, and, for gamma-distributed hit rates with this spread, a probability to miss a
lighter codeword at least 3 times larger than with equal hit rates. Then $e^{-\langle n\rangle}$ is optimistic, and a
larger `min_hits` is advisable; this happens, e.g., for the surface-code circuit `examples/surf_d5_H.mmx`.

**Probability to find a given codeword** (information sets; compare the QDistRnd manual, Sec. 3.2). A RW step reduces
$H$ (rank $r$) with a random column order; a codeword of weight $w$ is found if exactly one of its positions is among
the $s = n - r$ non-pivot columns ($s = \dim\ker H$). For uniformly random information sets, this happens with the
probability

$$
P_1(w) = \frac{w \binom{n-w}{s-1}}{\binom{n}{s}}
$$

per step ($P_1 = 0$ for $w > r + 1$), and the codeword is missed in $S$ steps with probability
$(1 - P_1)^S \approx e^{-S P_1}$. Among the weights $d_{\min} \le w < d_{\max}$, the weight with the smallest $P_1$ is
used, and the number of RW steps for a 1% miss probability is reported. Only the steps with uniform random column
permutations are counted (not those with localized windows; with `ksub>0`, the estimate is not available), which makes
the estimate conservative (for circuit DEMs, often very much so); on the other hand, as the QDistRnd manual points out,
the information sets of a sparse matrix need not be equally likely.

With `debug&1`, the binary prints $\langle n\rangle$ (`# codewords accumulated: ...`), the data for the
information-set estimate (`# RW information sets: n=..., rank(H)=..., steps=...`), and the estimate itself together
with the non-uniform hit count warning. The Python wrapper (`--verbose`) prints both estimates for a RW upper bound
which is not certified, with the RW steps and time needed for a 1% miss probability, and a model check which compares
the observed average hits of the codewords of weight $d_{\max}$ with $S \, P_1(d_{\max})$.

Example (`examples/surf_d5_H.mmx`, a circuit-level matrix of the distance-5 surface code; one thread for a reproducible
result):

```text
$ ./src/dist_m4ri method=1 finH=examples/surf_d5_H.mmx finL=examples/surf_d5_L.mmx steps=10000 seed=1 threads=1 debug=1
...
# codewords accumulated: total=100, min_w=5: cws=100, total_hits=208, hits min=1, max=33, avg=2.08, ...
# RW information sets: n=1958, rank(H)=120, steps=10000 (uniform permutations: 5000), ksub=0, ...
# Warning: non-uniform RW hit counts of the 100 codewords of weight 5: hits min=1, max=33, avg=2.08, ...
#   some codewords are found much more often than others of the same weight (relative spread of hit rates 1.89):
#   a lighter codeword may be missed with probability ~0.58 rather than exp(-1.70)=0.18; consider a larger min_hits
#   information-set estimate (uniform random information sets, n=1958, rank(H)=120): a codeword of weight 4
#   is found with probability 0.000845 per RW step, missed in 5000 uniform steps with probability 0.015 (1% ...
1 5 10000
```

Here, 10000 steps are not enough for `min_hits=5` ($\langle n\rangle = 2.08$). The hit counts of the 100 logical
operators of weight 5 range from 1 to 33, far from the Poisson distribution ($p = 4.3 \cdot 10^{-116}$), and the
probability to miss a lighter codeword is estimated as 0.58 instead of $e^{-\mu} = 0.18$, where $\mu = 1.70$ is the
Poisson mean fitted to the zero-truncated hit counts. The information-set estimate for the weight $w = 4$, the hardest
to find, gives the miss probability 0.015, and 1% after 10890 RW steps in total (half of them uniform, since
$n \ge 500$). In fact, the distance is 5: `method=2` with `wmax=4` finds no codeword in 0.1 s.

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

- **`outC=[file.nz]`**: Saves the discovered codewords of weight up to $w_{\min} + \text{dW}$ ($w_{\min}$: the minimum
  weight found) in standard **NZLIST** format (heavier codewords kept only for the `min_hits` statistic are not
  exported). With `method=2` or `3`, once the distance $d$ is known, the CC rounds $w = d, \dots, d + \text{dW}$
  enumerate all such codewords (unless the `timeout` is hit; with a stop target `dstop=U` below $d$, only the rounds up
  to $w = U + \text{dW}$ run); with `method=1`, only the codewords found by RW are exported:
  ```text
  %% NZLIST
  % generated by dist_m4ri
  <weight> <col_1> <col_2> ... <col_weight>
  ```
  *(Indices are 1-based).*
- **`finC=[file.nz]`**: Reads initial codewords from a file to initialize $d_{\max}$ and the codeword hash table.
- **`dW=[int]`**: When set (e.g. `dW=1`), preserves and exports codewords of weight up to $w \le d_{\min} + \text{dW}$.
- **`maxC=[int]`**: Limits collection to at most `maxC` unique codewords (of weight up to $w_{\min} + \text{dW}$).

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
./src/dist_m4ri (version 0.12.0): calculate distance of a classical or quantum CSS code
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
  dmax=[int]         Known upper bound on distance (CC up to w=dmax-1 only) (0)
  dstop=[int]        Stop once CC certifies dmin>=dstop (not an upper bound) (0)
  dexp=[int]         Expected distance for method=3 thread allocation (alias: dest) (0)

Search limits and stopping criteria:
  steps=[int]        Maximum RW decoding steps / information sets (100000)
  wmax=[int]         Maximum cluster weight to search in CC (0=until bound/timeout)
  wmin=[int]         Stop immediately if cw with weight <= wmin is found (1)
  min_hits=[int]     Stop RW when lowest-wt cws are hit min_hits times on average (5)
  timeout=[sec]      Execution timeout in seconds, 0 for infinite (60.0)

Multithreading and RW optimization:
  threads=[int]      Max worker threads to use (0: auto CPU count) (0)
  ksub=[int]         EXPERIMENTAL, do not use: subspace dimension sampled from ker(H) for RW (0)
  kwin=[int]         Localized column permutation window size W (0: auto/hybrid) (0)

Codeword collection:
  outC=[file]        Export found min-weight codewords (up to +dW) to file (.nz format)
  finC=[file]        Import initial codewords from file (.nz format)
  maxC=[int]         Maximum number of codewords to collect (0 for unlimited) (0)
  dW=[int]           Collect codewords up to the minimum weight found + dW (0)

Extra parameters (see --morehelp for details):
  smax=[int] (0)         Max syndrome weight for confinement profile (0 to disable)
  noscan=[0|1] (0)       Expert, method 2: CC at w=wmax only (no lower bound!)
  start=[list] (-1)      Expert: CC from listed columns only, e.g. start=0,48
  cbeg/cend=[int] (-1)   Expert: CC start column range, for split runs
  nothrottle=[0|1] (0)   Disable thread throttling (also --no-throttle)
  chunk_size=[int] (0)   RW batch chunk size (0: auto, alias: batch)
  win_mode=[0|1] (0)     Window mode: 0=Tanner BFS, 1=index proximity
  cov_cws=[int] (100)    Number of lowest-wt cws used for the min_hits statistic
  refresh=[int] (0)      RW steps between ker(H) basis refreshes (experimental ksub only; 5000 if ksub>0)
  seed=[int] (0)         RNG seed [0 for time(NULL)]
  debug=[int] (3)        Debug bitmask (0: silent, 1: summary, 2: progress, ...)

Help options:
  -h, --help         Display this help message (commonly used parameters)
  --morehelp         Display full help with all parameter descriptions
  --version          Display program version
```

### Full Parameter Listing (`dist_m4ri --morehelp`)

Use `./src/dist_m4ri --morehelp` to view exhaustive parameter explanations, including all `debug` bitmask flags (see
[Debug Output](#debug-output-debugint) below).

### Debug Output (`debug=[int]`)

All diagnostic output goes to `stderr`; `stdout` only has the result line `dmin dmax rw_steps`. The `debug` bitmap
selects the diagnostic output, where lower bits are more informative. The value is a decimal or hexadecimal (`0x...`)
integer. Multiple `debug` arguments are OR-combined, where the first one replaces the default `debug=3`: e.g.,
`debug=0` alone is silent, while `debug=4` and `debug=0 debug=4` both give only the periodic status.

- `0`: silent, except for errors, warnings on the validity of the result (expert options, the experimental `ksub`,
  invalid `finC` codewords), and the confinement profile with `smax>0`.
- `1` (default) **summary**: input matrices (sizes, numbers of nonzeros, maximum row and column weights), warnings,
  the reason why the run ended (`# stopped after ...s: ...`, e.g., `bounds coincide: dmin = dmax = 5`, `timeout=60s
  reached during CC round w=9 (37.5% of start columns claimed)`, or `RW convergence reached`), codeword and hit
  statistics, RW information sets, codeword export.
- `2` (default) **progress**: the run plan (method, threads, RW steps, CC weights, timeout, rng seed), finished CC
  rounds, new upper bounds, RW convergence.
- `4` **periodic status** after 1, 2, 4, ... seconds, then every 60 s: the bounds, the RW steps and their rate,
  $\langle n\rangle$, and the progress of the current CC round (or why no CC round runs).
- `8` **thread allocation and timing**: throttling, the planning and re-planning of CC rounds, waits for RW, and the
  predicted vs measured CC work of each round.
- `16` **code parameters**: $\mathrm{rank}\,H$, $\mathrm{rank}\,L$ (or $\mathrm{rank}\,G$), and $k$, by dense
  elimination at startup (slow for large matrices).
- `32` **codewords**: the support (1-based columns) of each new lightest codeword (RW, CC, `finC`).
- `64` **arguments**: the command-line arguments, the rng seed, and the files read.
- `128` **matrices** $H$, $G$, $L$ (for $n < 150$, unless bit `2048` is set).
- `256` **codeword dump**: all codewords in the hash table with their hit counts (at the end).
- `2048`: no size cutoff for matrices, codeword supports, and per-weight codeword statistics.
- `4096`: `dist_m4ri_old` only, CC recursion and confinement hash traces.

Bits `512`, `1024`, and `8192` to `32768` are reserved; bits `65536` and above are reserved for the Python wrapper
(ignored by the binary). E.g., `debug=7` adds the periodic status to the default output, and `debug=11` the thread
allocation and timing.

**Before version 0.11.0**, a different numbering was used: the old default bits `1` and `2` are now split into `1`,
`2`, and `8` (`debug=11` gives similar output); the old `4` (arguments) is now `64`; the old `8` (`dist_m4ri_old`: a
line every 1000 RW steps) is now `4`; the old `16` (new codewords) is now `32` (the upper bound messages: `2`); the old
`32` (matrices and codeword list) is now `128` and `256`; the old `64` and `128` (hash traces) are now `4096`.
Repeated `debug` arguments with different values were rejected.

### CLI Examples

```bash
# 1. Classical linear code using 8 threads in bracketing mode (method=3 is default); the number of RW steps varies
$ ./src/dist_m4ri finH=./examples/c204H.mmx dest=10 steps=100000 threads=8 debug=0
8 8 6154

# 2. Stim Detector Error Model (DEM) with timeout and codeword export (rw_steps=0: CC found a codeword of weight dmin)
$ ./src/dist_m4ri fdem=./examples/surf_d3.dem dexp=3 timeout=10 outC=cws.nz threads=4 debug=0
3 3 0

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
- `compute_quantum_distance(H, G=None, L=None, ...)`: Distance of one sector of a quantum CSS code, i.e., the minimum
  weight of a codeword $c$ with $Hc = 0$ and $Lc \neq 0$ (with $L$ constructed from $H$ and $G$ if not given).
- `compute_css_distance(Hx, Hz, Lx=None, Lz=None, ...)`: Distance $d = \min(d_X, d_Z)$ of a CSS quantum code (two runs
  of the binary, one for each sector); both `Hx` and `Hz` are required (otherwise `ValueError`; use
  `compute_quantum_distance()` for one sector). With `outC="cws.nz"`, the $X$- and $Z$-codewords are saved to
  `cws_X.nz` and `cws_Z.nz`; with `finC="cws.nz"`, the files `cws_X.nz` and `cws_Z.nz` are read if they exist, and
  otherwise the same file is given to both runs (codewords which are not valid in a sector are skipped). A known upper
  bound `dmax` on $d$ is not an upper bound of either sector: it is passed to both sector runs as the stop target
  `dstop` (and it is not stored in the cache), and the returned $d$ is at most `dmax`. The two sectors are cached as
  the quantum codes `(H=Hx, G=Hz)` ($d_Z$) and `(H=Hz, G=Hx)` ($d_X$), with `L=Lx` / `L=Lz` if given, i.e., with the
  same cache records as `compute_quantum_distance()`.
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
- `has_noise(circuit)` / `add_noise(circuit, p=0.001)`: Inspects a `stim.Circuit` for noise instructions, and adds
  uniform circuit-level noise to a noiseless circuit: `DEPOLARIZE2(p)` after two-qubit gates (except classically
  controlled ones, e.g., `CX rec[-1] 5`), `DEPOLARIZE1(p/10)` after single-qubit gates and on idle qubits in each
  `TICK`, and bit or phase flips with probability `p` after resets and before measurements (`compute_dem_distance()`
  does this for a noiseless circuit; for the distance, only which errors are possible matters, as long as `pmin=0`).
- `read_sparse_vectors(filepath)`: Parses NZLIST files into lists of 0-based integer support indices.
- Distance caching: `enable_distance_cache()`, `disable_distance_cache()`, `clear_distance_cache()`,
  `get_cached_distance(..., start=None)`, and the `cache_file` argument of the `compute_*_distance()` functions (a
  persistent JSON file with the version `"__version__"`; the CLI uses `tmp_dist_cache.json` in the working directory
  unless `--no-cache` or `cache=FILE` is given). A cache file written by a newer version is ignored (with a warning)
  and is not overwritten. The CSS records written before version 0.12.0 (keys `css:...`) are converted into the two
  sector records: the lower bounds are kept, but an upper bound only if a stored codeword confirms it, or if it is
  smaller than that of the other sector (earlier versions could store a `dmax` given to `compute_css_distance()` as the
  upper bound of both sectors). Earlier versions ignore a cache file written by version 0.12.0 (with a warning), and
  they may overwrite it with their own results: use separate cache files for different versions.
- Stop target `dstop` of all `compute_*_distance()` functions (CLI: `dstop=U`): the search ends once CC has certified
  `dmin >= dstop`, unless a lighter codeword is found (see [Bracketing Mode](#3-bracketing-mode-method3-default));
  the bounds are then `[dmin, w, rw_steps]` with `dmin >= dstop` and the weight `w` of the lightest codeword found
  (`0` if none). A cached lower bound `dmin >= dstop` answers the request without a calculation.
- Expert CC options of all `compute_*_distance()` functions: `start` (list of CC start columns, e.g. `start=[0, 48]`;
  separate cache record `<key>:start=a,b,c`) and `trust_start=True` (CLI: `--trust-start`; accept the start-list
  results as valid and copy them to the main cache record). The binary options `noscan`, `cbeg`, and `cend` are
  disabled in Python (ignored with a warning to `stderr`); see
  [Restricting the CC Search](#restricting-the-cc-search-expert-options).
- The experimental RW option `ksub` of the `compute_*_distance()` functions and the CLI should not be used: with
  `ksub>0` (and `method=1` or `3`), the same warning as from the binary is written to `stderr`.
- Optional solver backend: `solver="codedistance"` (uses the `codedistance` library if installed).
- Debug output: `debug=N` (default: 0) is a bitmap, where the bits `1` to `32768` are passed to the binary (see
  [Debug Output](#debug-output-debugint); its `stderr` is then printed), and the higher bits are used by the wrapper:
  `PY_DBG_COMMANDS` (`65536`: the `dist_m4ri` command lines and run times), `PY_DBG_CACHE` (`131072`: cache
  messages), `PY_DBG_CIRCUITS` (`262144`: Stim circuit to DEM conversion details), and `PY_DBG_KEEP_FILES`
  (`524288`: keep and list the temporary files). In the CLI, the value may be hexadecimal (e.g., `debug=0x10003`),
  multiple `debug` arguments are OR-combined, and `--verbose` implies the bits `1`, `2`, `65536`, and `131072`.

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

# 3. CSS Quantum Code ([[150, 32, 6]]): d = min(dX, dZ), with the bounds [dmin, dmax, rw_steps] of each sector
dist, d_x, d_z = dist_m4ri.compute_css_distance(
    Hx="examples/QX150.mtx",
    Hz="examples/QZ150.mtx",
    threads=4
)
print(f"CSS distance: {dist}, dX: {d_x}, dZ: {d_z}")  # 6, [6, 6, 0], [6, 6, 0]
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

# Run full C test suite (95 tests)
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


