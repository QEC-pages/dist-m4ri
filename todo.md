# Todo Notes for `dist_m4ri` Program

## 1. High-Priority Roadmap: High-Performance RW Overhaul

This work optimizes the Random Information-Set Decoding (RW, `method=1`) algorithm to achieve linear multicore scaling
(32–128 cores), minimal per-core memory footprint, and enhanced hit rates for quantum CSS codes and Detector Error
Models (DEMs).

Inspired by `sqetch` (arXiv:2607.28795, Appendix H) and the empirical convergence criteria in `QDistRnd`
(http://github.com/QEC-pages/QDistRnd).

### Implementation Strategy

- **Stage 1 (Immediate / Functional Baseline)**:
  - Implement Tasks 1–4 using `libm4ri` matrix primitives (`mzd_kernel_left_pluq`, `mzd_copy_row`, `mzd_echelonize`).
  - Validate mathematical correctness, hit-rate scaling, and backward compatibility (`ksub=0` retains original RW).
- **Stage 2 (Hardware-Optimized SIMD Engine)**:
  - For worker thread hot loops ($k_{\text{sub}} \le 64$), replace `libm4ri` with a zero-allocation, thread-safe,
    header-only AVX2 / AVX-512 bit-matrix kernel.
  - Utilize bit-sliced column representations (`uint64_t cols[n]`) to eliminate column transposition passes and
    completely remove the global `m4ri_mem_mutex`.

---

### Task 1: Subspace Sketching (`ksub`, experimental) for In-Cache Multicore Scaling
- [x] Add CLI parameter `ksub=[int]` (default `0` retains original full-matrix RW).
- [x] **Main Thread Precomputation**:
  - Compute a basis of the null space $N = \ker(H)$ (dimension $\nu \times n$, where $\nu = n - \mathrm{rank}(H)$) once
    using `mzd_kernel_left_pluq`.
  - Store $N$ as a global read-only matrix shared across worker threads.
- [x] **Worker Thread In-Cache Working Set**:
  - For $k_{\text{sub}} > 0$, allocate a compact $k_{\text{sub}} \times n$ matrix instead of the full $m \times n$
    check matrix.
  - For $k_{\text{sub}} = 64$ and $n = 5000$, working memory per thread is $\approx 40\text{--}80 \text{ KB}$, fitting
    entirely within private L1/L2 cache and eliminating L3 cache and DRAM bandwidth contention.
- [x] **Per-Trial Hot Loop**:
  1. Sample $k_{\text{sub}}$ distinct rows of $N$ (partial Fisher–Yates shuffle of the row indices; rows drawn with
     replacement would repeat and reduce the rank) into the thread-local working matrix `M_sub`.
  2. Perform in-place RREF via `gauss_one(M_sub, perm->values[i], rank)` in permuted column order (stopping as soon as
     `rank == ksub`, avoiding matrix transpositions and non-thread-safe `m4ri_mmc` allocations).
  3. Directly extract candidate codewords from the reduced non-zero rows (already in sorted original column
     coordinates; no dual transposition, coordinate unpermutation, or `qsort` needed).
  4. Test logical non-triviality against $L$ (column masks of $L$, see **Faster Logical Check**), and atomically
     update `dmax`.
- Findings (Oct 2026, `ksub` forced also for $m < \nu$, 8 threads): `ksub` explores only the span of $k_{\text{sub}}$
  rows of the fixed basis $N$ in each step, so a codeword is found only if all basis rows in its expansion are
  sampled. The hit counts are then concentrated on a few codewords, and `min_hits` stops RW early: e.g., for
  `bb72_si1000_X.dem`, after $\sim 2500$ steps with only 4–7 distinct codewords of weight 6 (full RW: 100 distinct
  codewords after $\sim 40000$ steps); for `torus60_d14_X.dem`, a codeword of weight 6 took 3–14 s instead of 3–15 ms.
  For `c1920H.mmx` (where `ksub` is allowed), the lightest codewords found in 20 s were similar to full RW. The
  fallback rule $m < \nu \Rightarrow$ `ksub=0` thus stays (it disables `ksub` for DEMs and most quantum codes).
- [x] `ksub` is kept for now, but marked experimental (should not be used) in the documentation, `--help`, and
  `--morehelp`; a warning is printed whenever `ksub>0` is given for RW (with or without `min_hits`, regardless of
  `debug`; in Python, once per `compute_*_distance()` call), as a reminder that the option does not do what it is
  supposed to.
- [ ] Find a subspace sampling for which the hit counts are uniform (or drop `ksub`).

---

### Task 2: Localized & Windowed Column Permutations for DEMs
- [x] Add parameters `kwin=[int]` (window size $W$, default `0` = uniform random permutation; alias `win=[int]`) and
  `win_mode=[0|1]` (default `0`).
- [x] **Window Construction**:
  - Pick a seed column $j_0 \in [0, n-1]$ uniformly at random per trial.
  - **Mode 0 (Tanner Graph BFS Neighbors)**: Perform BFS on the bipartite graph between columns and check rows starting
    from $j_0$ to gather $W$ geometrically/physically connected error mechanisms. Use a monotonic visit marker to avoid
    clearing arrays across trials.
  - **Mode 1 (Spacetime / Index Proximity)**: Define a contiguous coordinate window $[j_0, j_0 + W) \pmod n$ around
    $j_0$.
- [x] **Permutation Layout**:
  - Place the internally shuffled $W$ window columns in the first $W$ positions ($0 \dots W-1$).
  - Place the remaining $n - W$ columns at the tail ($W \dots n-1$).
- [x] **Support Overlap Filter & Windowed Echelonization**:
  - Optionally restrict row sampling to basis vectors of $N$ whose support intersects window $W$.
  - Exploit left-to-right echelonization: Gaussian elimination pivots first on the window columns, producing localized
    candidate codewords with support concentrated around the seed column $j_0$.

---

### Task 3: Adaptive Basis Evolution & Periodic Refresh of $N$
- [x] **Density Degradation Prevention**:
  - Avoid unconstrained random row mixing ($N_i \leftarrow N_i \oplus N_j$), which quickly drives basis row weights to
    $\approx n/2$ and degrades ISD hit rates.
- [x] **Steinitz Exchange with Discovered Codewords**:
  - Maintain an atomic pool of the lowest-weight non-trivial logical codewords found by worker threads.
  - Periodically substitute these lighter codewords into the basis $N$ using Gaussian elimination to maintain linear
    independence while systematically reducing the average row weight of the basis.
- [x] **Periodic Re-Echelonization (`refresh=[int]`)**:
  - Every `refresh` trials, generate a fresh canonical basis of $\ker(H)$ by re-echelonizing under a new global column
    permutation $\Pi$ (protected by a per-batch `pthread_rwlock_t`). This breaks sampling blind spots without
    accumulating row density.
- [ ] **Lock-Free Double Buffering (RCU Pattern)**:
  - Maintain `_Atomic(mzd_t *) N_active` and background buffer `N_staging`.
  - Worker threads acquire a local snapshot of `N_active` per batch without locking.
  - The coordinator prepares `N_staging` and atomically swaps `N_active` with release semantics.

---

### Task 4: Hit-Count Stopping Criterion (`QDistRnd` Convergence Style)
- [x] Add parameters `min_hits=[int]` (default `5`, `0` = disabled) and `cov_cws=[int]` (size of the representative
  set of lowest-weight codewords, default `100`).
- [x] **Codeword Multiplicity Tracking**:
  - Leverage the existing `cw_vec_t->cnt` field in `p->codewords` (which increments every time a duplicate codeword is
    found), with per-weight histograms of the numbers of codewords and of hits (`cw_cnt_w`, `cw_hits_w`).
- [x] **Representative Set**: keep up to `cov_cws` lowest-weight codewords (whole weight classes); heavier codewords
  are gradually replaced as lighter ones are found, regardless of `dW`; a new minimum weight is always recorded.
- [x] **Empirical Convergence Condition**:
  - In RW (`method=1` or `method=3`), stop when the average number of hits per codeword $\langle n\rangle$ reaches
    `min_hits`, both for the minimum-weight codewords and for the representative set (a single codeword suffices,
    e.g., for a classical code with $k=1$). As in QDistRnd, a lighter codeword is missed with probability
    $\approx e^{-\langle n\rangle}$, assuming that lighter codewords are found at least as often.
- [x] **Hit Uniformity Check**: at the end (`debug&1`, result not certified), a dispersion test of the hit counts in
  each weight class of the representative set (zero-truncated Poisson); a warning if they are strongly non-uniform.
- [x] **Information-Set Estimate** (Python `--verbose`, and in the binary with the hit uniformity warning): with
  uniform random information sets, a codeword of weight $w$ is found per RW step with probability
  $P_1(w) = w\binom{n-w}{s-1}/\binom{n}{s}$, $s = n - \mathrm{rank}\,H$; also the RW steps / time for a 1% miss
  probability (and the extrapolation of $e^{-\langle n\rangle}$ to 1%).
- [ ] **Proper estimate for matrices with many small-weight dual vectors**: see Task 6.

---

### Task 5: Native SIMD Bit-Matrix Engine (Stage 2)
- [ ] Implement a header-only, self-contained SIMD kernel (`src/simd_matrix.h`) for $k_{\text{sub}} \le 64$:
  - Flat bit-sliced column layout: `uint64_t cols[n]`.
  - Contiguous row layout: `uint64_t rows[64][(n + 63) / 64]`.
  - Permutation via single-instruction array indexing: `cols_perm[j] = cols[perm[j]]`.
  - Gaussian elimination via AVX2 (`_mm256_xor_si256`) and AVX-512 (`_mm512_xor_si512`) vector XOR loops.
  - Zero dynamic heap allocation per trial (`malloc`/`free` free inner loop).
  - 100% thread-safe; completely eliminate `m4ri_mem_mutex` contention in worker threads.

---

### Task 6: Miss-Probability Estimate for Matrices with Many Small-Weight Dual Vectors
- [ ] Find a proper estimate of the probability that RW misses a lighter codeword when the check matrix has many
  small-weight dual vectors (small-weight trivial vectors in $\ker H$). This is particularly bad for circuit DEMs,
  whose matrices have a huge number of column triplets summing to zero, as a consequence of the circuit structure: a
  weight-one error after a CX gate may give a weight-two error, and the net action of the three single-qubit errors
  is trivial. Both current estimates (Task 4) stay in place until then.
- Observations (Oct 2026, `method=1`, `min_hits=5`, `cov_cws=100`):
  - The uniform information-set estimate is accurate for some codes (c204: 166 expected vs 173 observed hits of the
    weight-8 codeword), but extremely pessimistic for DEMs: for `gross_uniform_X.dem` ($n=8784$, $\mathrm{rank}\,H =
    930$), $P_1 \approx 1.5\cdot 10^{-8}$ per step for $w=10$, i.e., $\sim 6\cdot 10^{8}$ steps for 1%, while the
    hit-based extrapolation needs $\sim 5\cdot 10^{3}$ steps. Steps with localized windows, which find local
    codewords much more often, are not counted in the uniform model.
  - The hit-based estimate with the representative set is conservative for codes with few minimum-weight codewords
    (c204: RW stops after $\sim 6000$ steps; the information-set estimate gives 1% after $\sim 90$ steps).
  - Hit counts can be strongly non-uniform within a weight class (`surf_d5`: some weight-5 logicals are hit 50–70
    times, others once), so that $e^{-\langle n\rangle}$ is optimistic there.
- Possible directions:
  - [x] Separate hit statistics for the uniform and the windowed RW steps (0.14.3): with `debug&1`, the binary prints
    the distinct codewords and hits of the three lowest weights found by each kind of step (`# RW step types, ...`),
    and the Python model check of the information-set estimate uses the hits of the uniform steps only. Codewords
    which violate the window model (e.g., time-like codewords with one error per round) are still found only by the
    uniform steps, which are kept for this.
  - Empirical per-weight hit rates from the representative set, extrapolated to weights below the minimum found.
  - Account for the small-weight trivial vectors explicitly (see **Trivial Codeword Statistics** below); this may
    also give a more accurate round estimate (degenerate multiplicity).
  - To check: the `ksub` basis refresh (Task 3) puts found low-weight codewords into the basis $N$, which may inflate
    their hit counts and $\langle n\rangle$.

---

## 2. Additional Algorithmic Explorations & Enhancements

- [x] **RW with pairs of non-pivot columns (Lee–Brickell / Stern, $p=2$)** (0.14.4, option `pairs`, on by default;
  `scratch/speedup/RESULTS.md`, experiment (h)): a codeword with exactly two non-pivot columns is the sum of their two
  circuits. With uniform information sets, a given codeword of weight $w$ is found
  $P_2/P_1 = (w-1)(n-r-1)/(2(r-w+2))$ times more often per step ($r = \mathrm{rank}\,H$, all pairs): about 10 for the
  dgross CSS code ($w = 18$), 20–100 for the reduced DEMs.
  - Prototype (numba, `scratch/speedup/p2_proto.py`, `p2_stern.py`): 10–35x more minimum-weight hits per step on
    bb_72, bb_90, and gross; all pairs take about 3x the time of the elimination, and the circuits are sparse, so that
    the buckets of a Stern collision search ($l$ random pivot rows) are very uneven.
  - In C: the non-pivot columns are sorted (radix sort) by a hash of their bits in $l$ random pivot rows, then by
    "zero L syndrome" and by weight. In a bucket, only the pairs whose weights differ by less than the weight limit are
    checked (word-wise popcount with early exit); the pairs of two columns with zero L syndromes (light trivial
    codewords, most columns of a circuit DEM) or with equal L syndromes are skipped (the rows of L reduced with the
    pivot rows, for at most 64 rows of L). The codewords of the pairs are recorded up to the weight of the upper bound
    (plus dW when collecting). Auto: $l = \lceil\log_2(n-r)\rceil$, reduced while a codeword of the weight of the upper
    bound is kept with the probability $\binom{r-l}{w-2}/\binom{r}{w-2} < 0.8$. The information-set estimate
    includes the $P_2$ term, and `debug&1` prints the pairs checked per step and the codewords from pairs and single
    columns.
  - Results (fixed RW steps): the minimum-weight codewords of circuit DEMs 7–26x more often per second (bb_72 7x,
    gross 22x, surf_d7 and mitten_150 26x); dgross: weight 18 in 3000 steps on 32 threads (20 without pairs; 2% more
    time per step); bb_360: 24 in 600 steps (26); classical LDPC codes with large distances: c1920H 210 instead of 320
    in 1000 steps (auto $l = 1$). Small codes: 2–5x longer steps (more hits recorded under the codeword mutex), but
    `min_hits` is reached in 2–34x fewer steps, in about the same time or less.
  - Possible follow-ups: thread-local batches of the recorded hits (the codeword mutex limits small codes with many
    threads, also without pairs); a cost model for $l$ (time of the pairs vs the elimination) instead of the fixed
    keep probability 0.8; triples ($p=3$) for small ranks.
- [x] **RW thread throttle for huge matrices** (0.14.3; `scratch/speedup/RESULTS.md`, experiments (a) and (e)): above
  15 MB of dense RW memory per thread, at most 32 RW threads, within a quarter of the physical memory (at least 2;
  1.5 GB if the physical memory is not known), instead of about 1.5 GB in total. The memory bus saturates at 32-64 RW
  threads whatever the matrix size.
  - bb_756 (reduced, 0.43 GB per thread): 32 RW threads and 1.49 steps/s by default, against 3 threads and 0.23
    steps/s before (6.5x).
  - bb_360 (49 MB per thread): 32 RW threads as before (45 steps/s).
- [x] **RW windows for big DEMs** (0.14.3; `scratch/speedup/RESULTS.md`, experiments (b), (e), (f)): the hybrid
  windows (every other RW step) have 256 columns for a DEM, 512 for matrices (`win_hyb`).
  - Windows of one and two rounds of columns (round estimate from the DEM) were slower than 512: a window of W
    consecutive DEM columns touches only so many detectors, e.g., one round ($L \approx 10D$ columns) about $2D$, so
    that most of its columns cannot be pivots.
  - A scan of single window sizes on six reduced DEMs (gross, dgross, bb_360, three mitten circuits): the best size
    is 192-256 columns for all of them (1.3-3.6 D), 2-10x better than 512; on dgross, about 10 distinct weight-18
    codewords per $2 \cdot 10^4$ steps with 256 columns, 0-2 with 512 or with one round. Cycling several sizes is
    worse than the best single size (the optimum is sharp). Matrices of codes gain nothing from 256.
  - The round estimate (D detectors and L columns per round) is printed with `debug&1` only.
  - Next:
    - the window size for much larger codes (bb_756: not measured; the optimum did not grow with D from 60 to 180);
    - an adaptive window size (e.g., from the light codewords found by each size), or a model of the optimum from
      the span of the minimum-weight codewords in the column order;
    - fewer uniform steps for big DEMs (keeping some for the statistics and for codewords which are not local);
    - window-restricted elimination (section 4).
- [ ] **Upper bounds from code logicals**: in a memory circuit, a minimum-weight logical of the code applied as
  single-qubit data errors at one time step flips no detector, so it is a DEM codeword of weight $d_{\rm code}$ (with
  si1000-type noise these errors are error mechanisms).
  - Python could map the data errors to DEM columns by their detector signatures and pass the codeword as `finC`.
  - This gives at once dgross $\le 18$, bb_360 $\le 24$, bb_756 $\le 34$ (RW: 18, 23, 58).
- [ ] **Incremental information sets (Canteaut–Chabaud)**, for huge DEMs: swap one column in or out of the
  information set per step (a single pivot change, about $r n/64$ word operations) instead of a full elimination.
  - This gives many more information sets per second (bb_756: 12–50 s per full elimination).
  - Consecutive sets are correlated: a separate upper-bound mode without `min_hits` and the miss estimate.
- [ ] **Trivial Codeword Statistics**:
  - Maintain counters of trivial codewords (stabilizers in quantum codes) of weights below $w_{\max}$ to monitor code
    degeneracy profile and ISD efficiency.
- [x] **Sorting Optimization**:
  - The codeword supports `ee` are sorted with the type-specific `rci_quick_sort` from `sort.h` (bitonic sort for up
    to 16 elements) instead of the generic `qsort`; with `ksub > 0`, the supports are extracted already sorted.
- [ ] **Sparse Gaussian Elimination**:
  - Evaluate whether sparse elimination (e.g. CSR row combining) can outperform dense bit-matrices for very large,
    highly sparse DEMs where $k_{\text{sub}}$ is small.
- [x] **Faster Logical Check**:
  - `sparse_syndrome_non_zero(L, ...)` cost $O(\mathrm{nnz}(L))$ per candidate, and it was the largest part of a RW
    step for DEMs, where most candidates below the weight limit are trivial (e.g., 248 candidates per step for
    `bb72_si1000_X.dem`, 1.5 ms of a 1.9 ms step), and also of CC rounds on DEMs (trivial zero-syndrome clusters).
  - Now the column bit masks of $L$ (`colmask_t`, `p->maskL`, $\lceil k/64\rceil$ words per column) give the check
    in $O(w)$ operations, for RW (also before sorting the support), `ksub` RW, and CC. Same results, one thread, the
    same RW steps (times include reading the input): `bb72_si1000_X.dem` RW 3.6 s $\to$ 0.63 s, CC (`wmax=6`)
    2.4 s $\to$ 1.1 s, `method=3` 3.2 s $\to$ 1.2 s; `gross_uniform_X.dem` 3.3 s $\to$ 1.5 s; `QX900` 0.68 s $\to$
    0.23 s; `surf_d5` 0.77 s $\to$ 0.51 s.
- [x] **CC row order** (`sort_rows`, 0.14.1): CC extends a cluster by the columns of the first unsatisfied row (in
  row order), so the binary sorts the rows of $H$ by increasing weight (stable) for matrices (`sort_rows=1`, the
  default; `sort_rows=0` for a DEM, whose detectors are in time order).
  - Extended CSS matrices with hook columns (`tmp/distance.py`; rows of weight 5–103): 1.2–2.3x to certify the
    distance, 1.7–4.2x to enumerate all codewords of a weight.
  - DEMs, and matrices exported from DEMs: 1.3–4x slower when sorted, since the first detector in time order has most
    of its columns before the start column of the cluster.
  - RW is unchanged (the same codewords with the same seed).
  - Choosing the row with the fewest admissible columns dynamically gave 2–23% fewer CC loop iterations on DEMs in a
    model, about as much as its per-node cost (`scratch/rowsort/`).
- [x] **RW Step Profile** (Oct 2026, one thread, after the faster logical check): the Gaussian elimination takes
  70–95% of a RW step for DEMs (`bb72_si1000_X.dem`: 215 of $\sim$310 µs; `gross_uniform_X.dem`: 6.8 of 7.8 ms), the
  transposition 7–15%, the extraction of the candidates the rest. Faster elimination (Task 5, or M4RI on a
  column-permuted copy, which needs thread-safe allocations) is the next step.
  - The weight of each candidate is now counted first (word-wise popcount, ending at the weight limit), and only the
    candidates below the limit are formed: `surf_d5` $-11\%$, `c1920H.mmx` $-9\%$ of the RW time.
  - With redundant rows of $H$, the RW matrices keep only the $\mathrm{rank}\,H$ non-zero rows after the first step,
    so that the transposition skips the zero rows.
  - Contention on `cw_mutex` is negligible: only the candidates below the weight limit and with $Lc \neq 0$ are
    recorded, and with `min_hits` about `min_hits` $\times$ `cov_cws` hits are recorded until convergence (e.g., 1.4 µs
    of lock wait per step of $\sim$260 µs for `surf_d5` with 32 threads); it only shows in runs of a few
    milliseconds.
- [x] **Distance Benchmark Update**: the overnight run of 0.14.0 (Oct 2026, `benchmark/BENCHMARK.md`) certified the
  distances of all benchmark circuits up to gross. Still open: double gross, bb_360, and bb_756; see **Bounds for the
  largest circuits** in section 4.

---

## 3. Operational Modes Reference

### Mode: Distance Verification
- Given a suspected distance $d$ (e.g., $d = 10$), run `method=1 dmax=10 wmin=9 min_hits=0`.
- RW ignores codewords of weight $\ge 10$ and ends as soon as it finds a non-trivial codeword of weight $w \le 9$
  (output `1 w S`; with `debug&1`, the stop line `found a codeword of weight w <= wmin=9`).
- If no lighter codeword is found in the `steps` RW steps, the output is `1 10 S` (the supplied upper bound).
  (The legacy `dist_m4ri_old method=1 wmax=10 wmin=9` outputs `-w` or `0`.)

### Mode: Upper Bound Search with Hit-Count Saturation
- Run RW with `method=1 steps=N min_hits=M cov_cws=K`.
- Collects the lowest-weight codewords into the hash table, terminating when the average number of hits per codeword
  reaches $M$, both for the minimum-weight codewords and for the representative set of $K$ lowest-weight codewords.

### More debug bits
- [x] Debug bits reorganized by usability (see README "Debug Output" and `--morehelp`): 1 summary (input summary,
  stop reason), 2 progress (run plan with the seed), 4 periodic status (1, 2, 4, ... s, then every 60 s), 8 thread
  allocation and timing, 16 code parameters, 32 codeword supports, 64 arguments, 128 matrices, 256 codeword dump,
  2048 no size cutoff, 4096 legacy traces; bits from 65536 for the Python wrapper (command lines, cache, circuits,
  keep temporary files). Repeated `debug=` arguments are OR-combined (hexadecimal values accepted).
- [x] Detailed timing information (measured / predicted values, reasons for termination): debug bits 8 and 1.
- [x] ranks of the matrices and their sizes and code dimensions (to avoid doing rank in python)
  - [x] After RW, with `debug&1`, the binary prints `# RW information sets: n=..., rank(H)=...,
    steps=... (uniform permutations: ...), ksub=..., ... s/step per thread, ... RW threads` (parsed by Python).
  - [x] With `debug&16`, the binary prints rank(H), rank(L) or rank(G), and k (dense elimination at startup); sizes,
    numbers of nonzeros, and maximum row and column weights are in the input summary (`debug&1`).
  - [ ] Python could parse the `# code parameters:` line instead of computing ranks itself. Needed for
    `tmp/distance.py`, e.g., an option of the binary to report the ranks, $k$, and other parameters to Python.
- [x] in addition, print estimated prob to find (miss) a single codeword of given
  weight based on the number of information sets and the n and distance found with RW.
  - Done: Python `--verbose`, and the binary with the hit uniformity warning (Task 4); a proper estimate for DEMs
    is still open (Task 6).
- [ ] The periodic status shows that in `method=3` a CC round can wait for a whole RW batch before any CC thread
  starts (e.g., `c1920H.mmx`: 100 RW steps of 6.6 ms per batch, i.e., 0.66 s): smaller RW batches while CC rounds
  are short, or an RW batch interruption when a CC round starts.
  - The batch size is set in steps (25–50 for $n \ge 5000$), not in time. At 0.7 s per step (bb_360) a batch takes
    20–35 s, and at 12–50 s per step (bb_756) several minutes.
  - Threads which cannot run RW (beyond the RW throttle) join a CC round at once, so the wait matters mainly when
    all threads run RW.
  - Proposal: batches of at most about 0.25 s by the measured step time, or leave the batch when a CC round needs the
    thread.

### More improvements 
- [x] try verifying (using cluster decomposition from vecdec (?)) if the
      constructed / imported non-trivial CW can be decomposed into disconnected
      components.  Keep only non-trivial part (does it ever happen? how
      expensive is this?)
  - Answer (Oct 2026): it never happens for RW. A RW candidate is the fundamental circuit of a non-pivot column (the
    unique kernel vector on the pivots plus that column), i.e., a minimal codeword, and a disconnected codeword is never
    minimal (each component is itself in $\ker H$). Instrumented test: 0 disconnected among $\sim 7\cdot 10^6$ kernel
    vectors (surf_d5, bb72 and dgross DEMs, mitten_90 and mitten_270 full DEMs, dgross CSS, c1920, a random LDPC code).
    Only the experimental `ksub` rows are not circuits; CC clusters are connected by construction.
- [ ] See that matrices constructed when generating L from H (e.g.) are reused when codeword searching.
- [ ] For a very large circuit, when RW is very slow, should we try constructing
      codewords from the underlying code (can this even be done from a DEM ---
      or should it be done in python?)
  - See **Full DEM from the simplified DEM** below (done by hand for mitten_270 in Oct 2026).
- [ ] 

---

## 4. Large DEMs, Simplified DEMs, and Symmetry (discussion, Oct 2026)

- [x] **No equal columns in DEMs**: stim's DEM of a circuit with a REPEAT block has duplicate error mechanisms even
  after flattening the DEM (loop folding), e.g., mitten_90 custom: 23394 errors, 22170 distinct (the DEM of the
  flattened circuit has none); about 16% of the columns of the BB DEMs. CC treats each duplicate separately: merging
  them gave 1.6–1.7x for CC (dem150 and dgross DEMs, `wmax=6`, 128 threads), and RW eliminates fewer columns.
  - Python: merge equal columns (XOR-reduced detector and observable sets, combined probability $p_1 + p_2 - 2p_1p_2$)
    when writing the DEM (or flatten the circuit before generating the DEM; this loses stim's loop folding for very
    long circuits). Not needed now that the binary merges (and it would change the DEM cache keys).
  - C: done, the binary merges the equal columns on reading a DEM (`merge=1`, the default), with a map back to the
    original columns (`outC`, `finC`, `start`, `cbeg`, `cend`); bb72 DEM: 2664 -> 2232 columns, CC `wmax=6` in 0.06 s
    instead of 0.21 s.
- [x] **Full DEM from the simplified DEM** (`--full` mode of CSS circuits): the simplified DEM is the full DEM with the
  minority-basis detector rows deleted (equal columns merged, empty columns dropped; they never flip an observable), so
  $d_{\rm simp} \le d_{\rm full}$, and $d_{\rm full} = d_{\rm simp}$ iff some minimum-weight simplified codeword lifts
  (one full error per column, with the deleted detectors cancelling too). Sufficient and checked column by column:
  every simplified error has a "clean" preimage (no deleted detector), which holds for CSS circuits with depolarizing
  noise. Checked for mitten_90 and mitten_270 (custom, hef): all simplified errors have clean preimages (6.4–6.7
  preimages each); RW on the simplified DEMs found weight-7 (hef) and weight-9 (custom) codewords in about 2 s, which
  lift to the full DEMs, so mitten_270 full: hef $d = 7$, custom $d = 9$ (RW on the full DEMs had not found them in
  2.8 h: the per-step probability scales like $(\mathrm{rank}/n)^{w-1}$).
  - Python: build both DEMs with a detector map (`strip_minority_detectors` returning the kept detector indices),
    check the clean preimages, run (or take from the cache) the simplified result, lift its codewords to the full DEM
    (`finC`, upper bound), and use $d_{\rm simp}$ as `dmin` for the full run; run the full binary only if a gap remains.
  - Default on for `--full`, but only for circuits verified as CSS: the circuit classification, plus the DEM-level
    checks (no full error projects to an observable-only column; every simplified error used has a clean preimage),
    and every lifted codeword re-verified in the full DEM (`finC` validation). The lower bound
    $d_{\rm simp} \le d_{\rm full}$ alone holds for any set of deleted detectors. An option disables it (for tests).
  - Done differently (exact reduction in the binary, `drop=FILE`): every full column touching a deleted detector needs
    a twin (same kept detectors, no deleted one, observables $L + A H$, $A$ by elimination); then the binary runs on
    the kept detectors and the columns touching no deleted detector, and its codewords are full codewords (no lifting).
    Python passes the minority-basis detectors with `--full` (`drop=1`, the default). mitten_270 full: 70254 -> 8910
    columns, RW finds weight 9 (custom) and 7 (hef) in 3.5 s; mitten_90 full: exact $d = 5, 4$ by CC in 0.1 s.
  - Since 0.13.0, the Python options `--simple` and `--full` (`simple=`, `full=`) are disabled (ignored with a
    warning): the DEM of a Stim circuit always has all detectors.
  - Done in 0.14.0: `drop=auto` in the binary (the default, also `drop=1`; `drop=0` keeps all detectors) finds the
    largest set of detectors which can be dropped with $A = 0$ from the DEM alone (sub-columns with the same
    observables, groups of the O-minimal columns, seeds with propagation, see the README), and verifies it as for
    `drop=FILE`; it finds exactly the minority-basis detectors for all CSS circuits tested (bb_756: 12474 of 25704
    detectors in 1.3 s). `outD=FILE` writes the detectors dropped; Python: `drop=css` (the classification of the
    circuit; for a DEM, `drop=auto` with a warning) and `--out-reduced [PREFIX]` (the reduced DEM `PREFIX_red.dem` and
    its maps `PREFIX_red.json`).
  - Open: a second pass allowing $A \neq 0$ (observables which include detectors of the dropped set).
- [x] **Uniform noise model** (0.14.0): `add_noise()` classifies the gates with `stim.gate_data()` (all one- and
  two-qubit gates get depolarizing noise, a classically controlled gate `DEPOLARIZE1(p/10)` on its qubit; no noise for
  `MPP`, `MXX`, `MYY`, `MZZ`, `SPP`); `strip_noise()` (stim's `without_noise()`); `has_noise()` also detects noisy
  measurements such as `M(0.001)`; CLI `--p-noise P` (as in `tmp/distance.py`; 0: none added) and `--force-noise`
  (replace the noise of a circuit, e.g., to compare circuits with the same noise model; with `--p-noise 0`, remove it).
- [ ] **Extended CSS distance as a seed** (`tmp/distance.py`, code-capacity matrices with one column per hook error):
  not guaranteed to relate to the DEM distances for every CSS circuit. Its lower-bound direction
  ($d_{\rm ext} \le d_{\rm simp}$, via the telescoping sum of the detectors of each check) needs every single fault's
  net data effect to be a data qubit, an included hook, or nothing, and its codewords lift only if the hook errors fit
  the timing of the other checks. Observed: custom $d_{\rm ext} = 4, 8 < d = 5, 9$ (mitten_90, mitten_270; the
  extended codewords cannot lift); hef $d_{\rm ext} = 6, 10 > d = 4, 7$ (the hook extraction finds no X-check
  ancillas in the hef circuits, `Hx` with 0 rows, so no hooks are added). A rigorous variant: build the model from
  the DEM itself (net data effect of each fault, e.g., with one extra observable per final data measurement), and
  check rank($[H_K; H_z N]$) = rank($H_K$) and $L = L_z N$; its distance is then a certified lower bound on
  $d_{\rm simp}$, and it is code-sized (and keeps the symmetry of a symmetric schedule).
  - Prototype of the rigorous variant (Oct 2026, `scratch/speedup/model_lb.py`).
    - Single data measurements are not deterministic, so the circuit gets one deterministic observable $O_h$ per
      check $h$ of the measured basis: the parity of the final data measurements on its support, read from the
      detectors after the final measurement.
    - Verify rank($[H; O]$) = rank($H$), i.e., the rows $O$ lie in the row space of $H$.
    - The model is $(O, L)$ with equal columns merged: a row reduction with $L$ unchanged, so its distance is a
      certified lower bound on the circuit distance.
    - Model sizes and distances (the verification held for all of them):
      - bb_72: 36 x 180, distance 5 ($d = 6$);
      - bb_90: 45 x 225, 7 ($d = 8$);
      - bb_108: 54 x 270, 6 ($d = 8$);
      - gross: 72 x 360, 8 ($d = 10$);
      - dgross: 144 x 720, distance 12 (no model codeword of weight $\le 11$, 550 s on 32 threads; one of weight 12
        found in 467 s on 96 threads), so $d \ge 12$ (the circuit DEM took 2.8 h on 128 threads to certify 11);
      - bb_360: 180 x 900, $\ge 12$ (no model codeword of weight $\le 11$, 966 s on 32 threads).
    - Not tight for these schedules (1–2 below $d$), but code-sized; the lifting and refinement below could close the
      gap.
- [ ] **Row reductions of a DEM, with lifting and verification in C** (generalizes the two items above): replacing $H$
  by $RH$ with rows in the row space of $H$ (deleting or summing detector rows), with the observables $L$ unchanged,
  can only lower the distance, $d_R \le d_{\rm full}$; deleting or merging observables could raise it. Columns equal
  on $(RH, L)$ merge into classes, and $d_R = d_{\rm full}$ iff some minimum-weight reduced codeword lifts (one full
  column per class, with $Hx = 0$).
  - Exact reductions ("clean twins"): a detector set $M$ can be deleted together with all columns touching it if each
    such column has a twin with the same pattern on the other detectors and on $L$ that does not touch $M$ (columns
    with no other detector and no observable are dropped). The reduced DEM is a submatrix of the full one with the
    same distance, and its codewords are full codewords (no lifting). The simplified DEM is the CSS case ($M$ = the
    minority detectors). The condition must hold for $M$ as a whole (single minority detectors usually fail it), so
    candidate sets (sectors, coordinate classes) are tested as groups.
  - Lossy reductions: summing the detectors of each check over blocks of rounds (over all rounds: the code-sized model
    of the net data effects, see above) gives cheap certified lower bounds, tight or not depending on the schedule.
    Lifting: one member per class with the deleted or summed rows cancelling (meet-in-the-middle over the members,
    within a time window), else a repair with extra columns (a heavier upper bound); refinement: split the summed rows
    of the checks where a lift failed back into per-round detectors and solve again (row generation), until a lift
    matches the lower bound. Merge the equal columns first: summed rows otherwise get huge row weights for CC.
  - In C: the full DEM and a reduction spec (a delete list or an aggregation map from Python, or rules on detector
    coordinates); build $RH$, merge equal columns (map from each class to its full columns), run RW or CC on the
    reduced matrix, lift and verify each candidate ($Hx = 0$, $Lx \ne 0$ for the full DEM), and write the full
    codewords to `outC`. Report $d_R$ (certified once CC is done) as `dmin` and the verified lifts as `dmax` for the
    full DEM; `min_hits` then refers to $d_R$, so exactness only when the bounds meet.
  - Order: clean-twin deletion together with equal-column merging (replaces the Python simplified/full detector map,
    verified automatically); then the code-sized model as a lower bound for large distances; adaptive lifting and
    refinement as research.
- [ ] **Bounds for the largest circuits** (X sector, Oct 2026, `scratch/speedup/RESULTS.md`):
  - double gross: $[12, 18]$, before $[11, 19]$. The lower bound comes from the code-sized model above. The upper bound
    comes from RW with index windows of 1 round (3 weight-18 codewords in 180 s, in `scratch/speedup/cws/`); they are
    errors of the final data measurements, i.e., a code logical at the time boundary.
  - bb_360: $[12, 23]$, before $[4, 25]$. The lower bound comes from the model; the upper bound from RW, the default
    hybrid on 16 threads, found twice and saved in `scratch/speedup/cws/`.
  - bb_756: only RW, lightest 58 (the code has distance at most 34).
  - With the model 1–2 below $d$ for the smaller BB circuits, the true distances may be well below the upper bounds
    (e.g., 13–14 for double gross if the gap is the same). Lifting and refinement of the model codewords, or RW with
    $p=2$ and better windows, could close the gap.
- [ ] **Symmetry in code mode** (not for DEMs: their observables break the symmetry, |Aut| = 1 for the BB DEMs once
  duplicates are merged). Automorphism groups with BLISS (python-igraph) take milliseconds (gross, double gross, a
  generalized bicycle code with $n = 3000$, QX900) to 0.26 s (the affine-plane incidence code AG(2,64), $n = 4096$,
  $|\mathrm{Aut}| = 4\cdot 10^{11}$), so they could be computed by default.
  - Use: one start column per column orbit (`start` list, unlimited clusters), when the orbit size is large compared
    with $w$ (gain about orbit size / $w$; measured 7.2x for double gross with `start=0,144`, one column per orbit of
    the translations; the full group, $|\mathrm{Aut}| = 864$, has a single column orbit). Orbit exclusion: the k-th
    start column only needs clusters avoiding the orbits of the previous ones (about 2x more for two orbits).
  - Correctness: the graph must contain all defining matrices ((Hx, Hz) or (H, G): an automorphism then preserves both
    row spaces); with `finL`, verify each generator by rank checks of $[H; L]$. Certify `dmin` without `trust_start`
    (e.g., the binary verifies the given generators and computes the orbits).
  - Codeword lists (`outC`, `maxC`, the Python cache): by default the shorter list of canonical orbit representatives
    (e.g., the lexicographically smallest image; the group recorded with the list, so that full lists from runs
    without symmetry can be canonicalized when merged); the expanded list only on request (it can be huge: e.g.,
    $|\mathrm{Aut}| = 4\cdot 10^{11}$ for AG(2,64)).
  - Dependency: python-igraph is optional (networkx is too slow for large graphs); or nauty/Traces in C.
- [ ] **Orbit-canonical RW hit counts** (with a known group): with uniform permutations, all codewords of an orbit have
  the same hit rate, and a lighter orbit $O'$ is missed with probability $\approx e^{-|O'| p N}$, with $|O'| \ge |G|/w$
  for a group acting freely on the columns. Count hits per orbit (uniform steps only), canonical form: cheap for
  translation groups (QC, BB, GB codes), group enumeration for $|G| \lesssim 10^5$, canonical labeling otherwise.
- [ ] **Time-translation classes of DEM codewords** (many rounds): group the translates of a codeword (detector
  coordinates or the round of each detector, kept when unrolling) for `outC` and the reported hit counts. Not for the
  miss probability: the symmetry is approximate, and a lighter codeword may sit at a time boundary (no translates).
- [ ] **RW window steps restricted to the window**: the circuits of the non-pivot window columns are those of the window
  submatrix (all rows touched by the window), at a fraction of the cost of the full elimination; window shapes: time
  slabs (space-like logicals), tubes of a few stabilizers over all rounds (time-like logicals), Tanner-graph balls, and
  the full matrix (uniform steps, which keep every codeword reachable). Separate hit statistics for the window steps.

