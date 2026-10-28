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
  - [x] 0.14.6: a work bound for the auto $l$ (`benchmark/BENCHMARK.md`, Section 8; `scratch/bench_0.14.5/`): $l$ is
    at least the $l$ at which the about $(n-r)^2/2^{l+1}$ pairs of a step are within about the work of its
    elimination, $rn/16$ pairs (or $2^{20}$ pairs for small matrices), and the stop flags and the upper bound are also
    polled within the buckets (every 65536 pairs). Early in a run, while the upper bound is far above the distance,
    the keep probability alone gave $l = 1$ (all pairs in two buckets): 0.14.5 on `bb_756` completed no RW step in
    the first 32 s and did not stop at its timeout, and on `mitten_270` took 0.6 s of RW to the distance instead of
    0.16 s (0.14.3). The bound gives $l \ge 6$ for these DEMs (0.14.6: 0.12 s of RW to the distance of `mitten_270`,
    and `bb_756` stops at its timeout), and no change for codes with $n - r < r$ (`c1920H`).
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
  - [ ] **Code distance of a circuit as a fast upper bound** (for DEMs and Stim circuits, especially large ones):
    extract the code from the circuit with the functions of `tmp/distance.py` (`extract_matrices_and_hooks()`: `Hx`,
    `Hz`, `Lx`, `Lz` from the ancilla interactions and the observables; `compute_css_distance()` for the code), and
    compute the distance of the sector of the memory basis (Z-type logicals for an X-basis memory: `H = Hx`,
    `L = Lx`). It is code-sized (seconds where the DEM takes hours) and an upper bound on the circuit distance whenever
    the logical can be applied by single faults (data errors, or flips of the final data measurements, as with the
    uniform noise).
    - Certify it: lift the logical to a DEM codeword (the error mechanisms of the flips of the final data measurements
      on its support, found by their detector signatures), check $Hx = 0$, $Lx \ne 0$, and pass it to the DEM run as
      `finC` (an upper bound `dmax` for CC and RW); report a note if no lift exists.
    - A DEM-only variant without the qubit classification of `distance.py` (which fails, e.g., for the hef circuits,
      with no X-check ancillas found, and for routing circuits): the error mechanisms which touch only the final
      detectors (the flips of the final data measurements) form a sub-DEM which is the code itself (rows: the checks
      of the measured basis, as in the code-sized model); its distance is an upper bound, and its codewords are DEM
      codewords.
    - Python: a function (e.g., `compute_code_upper_bound(circuit)`) and a CLI option, the bound stored as `dmax` in
      the cache record of the DEM (with its codewords); together with `compute_model_lower_bound()`, the bracket
      $[d_{\rm model}, d_{\rm code}]$ for the circuit without CC or RW on the DEM.
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
- [ ] **References for the acceleration techniques** (search online): find the prior work and the right citations for
  each technique used by `dist_m4ri` (for the README, `introduction.md`, and a paper), and whether better variants
  exist. Starting points to verify:
  - RW and information sets: Prange (1962); Lee–Brickell (1988) and Stern (1989) for the pairs ($p = 2$) and the
    collision buckets; Canteaut–Chabaud (incremental information sets); later information-set decoding (Dumer;
    May–Meurer–Thomae; Becker–Joux–May–Meurer) for further gains; QDistRnd (random information sets for quantum
    codes, the hit-count convergence criterion); prior use of localized windows or time slabs for circuit-level
    distances.
  - Elimination: the Method of the Four Russians (Arlazarov, Dinic, Kronrod, Faradzev) and M4RI (Albrecht, Bard;
    Gray-code tables), for the scattered column orders of RW.
  - CC: connected-cluster enumeration (Dumer, Kovalev, Pryadko 2017), the meet-in-the-middle hash for the last two
    levels, work units and thread scheduling; the branch-and-bound pruning (forbidden set, fewest-available row,
    packing bound) of J. Y. Lee, Okada, Maskara, Kasai, Zhou, arXiv:2609.30069, App. B, used by G. Lee,
    arXiv:2610.13880, App. C.4 (see **Branch-and-bound CC** below).
  - Statistics: index-of-dispersion tests and zero-truncated Poisson fits of the hit counts.
  - Circuit-level distance: stim's search for undetectable logical errors, SAT/MaxSAT and ILP approaches (e.g., for
    the BB codes), spacetime codes (Delfosse, Paetznick), hook errors and distance-preserving schedules, extended CSS
    distances with hooks, exact reductions of DEMs (dropped detectors), and code-sized or folded models as certified
    lower bounds.
  - Memory-bandwidth-aware thread throttling for dense GF(2) linear algebra.

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
  - [x] Implemented as `compute_model_lower_bound()` (0.14.5; CLI `--model-lb`, `--out-model PREFIX`): the rows
    are the sums of the chains of detectors linked by a shared measurement record other than a final data
    measurement, found from the circuit (one chain per check of a memory circuit), so that the rows lie in the row
    space of $H$ by construction (no rank check); the chains of the minority-basis checks are extra rows
    (`minority=True`). The binary runs on the measured-basis rows first, then RW and, if needed, CC on all rows from
    that bound; the bound is stored as `dmin` (and `"dmin_model"`) in the cache record of the DEM.
    - Results (16 threads, 120 s per run): tight for the stim surface-code circuits ($d = 3, 5$, X and Z), mitten_90
      hef (4), mitten_270 hef (7), and the torus n120 circuit (6); below $d$ for mitten_90 custom 4 ($d = 5$) and
      mitten_270 custom 8 ($d = 9$; the same as $d_{\rm ext}$), mitten_150 custom 6 ($d = 7$), the torus n60
      circuit 4 ($d = 6$), and bb_72 5 (X and Z), bb_90 7, bb_108 6, gross 8, as the prototype. In the torus circuits
      all qubits are measured in each round (the final data measurements are then all those of the final block).
    - The minority-basis rows raised none of these bounds within the time. CC on the model with all rows alone is
      much slower (bb_90: 53 s instead of 0.1 s; gross: only 7 certified in 75 s), which the RW stage avoids (gross:
      1.8 s in total). For mitten_270 custom, RW on all rows found no codeword of weight 8 (lightest 9, the circuit
      distance) and CC did not finish $w = 8$ in 120 s: a longer run may show whether the minority rows raise the
      bound to 9 there.
    - [x] 0.14.7: the minority-basis rows are off by default (`minority=True`, CLI `--model-minority`, to use them).
      They raised none of the bounds of the 16 circuits of `benchmark/BENCHMARK.md` Section 8.5 (`mitten_270` custom
      stays 8, with CC on all rows in 7.8 s) and of the stim surface codes ($d = 5, 7$), and cost up to 9 s: with
      depolarizing-type noise, every net data error which flips the observable of the measured basis also occurs
      without a component of the other basis, so that a lightest model codeword has a version which these rows do not
      detect. They may raise the bound for other noise models (e.g., biased noise).
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
  - Lossy reductions: summing the detectors of each check over all rounds gives the code-sized model of the net data
    effects (see above), a cheap certified lower bound, tight or not depending on the schedule. Sums over blocks of
    at least 2 rounds give exactly the same bound: a model codeword realized within one block (its partial syndromes
    at rounds $t$ and $t+1$ in the same row) is invisible again; see the folded models below.
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
- [ ] **Folded models (research)**: a refinement of the code-sized model which keeps the timing within the rounds.
  - Rows $R_{h,c} = \sum_{t \equiv c \pmod p} D_h^t$ for each check $h$ and residue $c$, with $p \ge 2$ (the boundary
    detectors in their residue classes): sums of detector rows, so that the bound is valid (any row reduction is).
  - Consecutive rounds are in different rows: the partial syndrome $P$ which the faults of one round leave (detectors
    $P$ at rounds $t$ and $t+1$, which cancel in the code-sized model) is seen as in the DEM, which is why the lightest
    codewords of the code-sized model do not lift.
  - Bounds: $d_{\rm model} \le d_{\rm fold}(p) \le d$, and $d_{\rm fold}(p) \le d_{\rm fold}(p')$ if $p$ divides
    $p'$ (each mod-$p$ row is a sum of mod-$p'$ rows).
  - Size: like a DEM with $p$ rounds and periodic time, e.g., about $300 \times 3000$ for double gross with $p = 2$
    (reduced DEM: $2736 \times 26208$). Spurious codewords remain (equal partial syndromes $p$ rounds apart cancel),
    and the row weights are those of the DEM, so that CC costs more than on the code-sized model.
  - Plan: the round of each detector from its position in its chain (`_model_detector_groups()`); test $p = 2, 3, 4$
    on bb_72 and gross (code-sized bounds 5 and 8, $d = 6$ and 10), then double gross and bb_360 if CC reaches the
    needed weights.
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
    start column only needs clusters avoiding the orbits of the previous ones (about 2x more for two orbits; details
    in **Orbit exclusion** below).
  - Correctness: the graph must contain all defining matrices ((Hx, Hz) or (H, G): an automorphism then preserves both
    row spaces); with `finL`, verify each generator by rank checks of $[H; L]$. Certify `dmin` without `trust_start`
    (e.g., the binary verifies the given generators and computes the orbits).
  - Codeword lists (`outC`, `maxC`, the Python cache): by default the shorter list of canonical orbit representatives
    (e.g., the lexicographically smallest image; the group recorded with the list, so that full lists from runs
    without symmetry can be canonicalized when merged); the expanded list only on request (it can be huge: e.g.,
    $|\mathrm{Aut}| = 4\cdot 10^{11}$ for AG(2,64)).
  - Dependency: python-igraph is optional (networkx is too slow for large graphs); or nauty/Traces in C.
- [ ] **Orbit exclusion** (for the programs which construct codes with known symmetries; needs the forbidden set of
  the **Branch-and-bound CC** below): the CC search from the $j$-th start column forbids all columns of the orbits of
  the earlier start columns.
  - **Partly done in 0.15.1:** `start_scan=1` (expert, with `start`): a cluster grown from a listed column $c$ only
    uses the columns $> c$, which is valid if the listed columns are the smallest columns of the orbits (one per orbit:
    every column below a listed one is in an earlier orbit), and it is exactly the orbit exclusion if each orbit is a
    range of consecutive columns, e.g., `start=0,l start_scan=1` for the two blocks of a 2BGA code over an abelian
    group (2x less CC work than `start=0,l`), or after renumbering the columns by orbits (largest first) in the program
    which builds the code. `scratch/bench_0.15.1_scan/`. Still open: the inputs below (orbits in any column order,
    verified generators), and the work-unit split for few start columns (`scratch/symmetry/REPORT.md`).
  - Input, one of:
    - `orbits=FILE`: an orbit label for each column ($n$ integers, in the original 0-based column order), and one start
      column per orbit (from `start=`, by default the smallest column of each orbit). Trusted, as `start`: a
      `# WARNING:`, and `dmin` is valid only if the labels are the orbits of a group of code automorphisms;
    - `aut=FILE`: generators of the group, each a permutation of the columns (a line of $n$ 0-based images). The binary
      verifies each generator $\pi$ by rank tests ($\mathrm{rank}\,[H; H\pi] = \mathrm{rank}\,H$, the same for $G$, or
      for $[H; L]$ with `finL`), computes the orbits (union-find over the cycles of the generators) and the start
      columns, and the result is certified without `trust_start`.
  - Search: start columns $r_1, \dots, r_k$ of the orbits $O_1, \dots, O_k$; from $r_j$, unlimited clusters (as with
    `start`) with $F = O_1 \cup \dots \cup O_{j-1}$ at the root. A codeword $C$ is found at most once: from the first
    orbit $O_j$ it meets, if $r_j \in C$. Every codeword is the image of a found one: $g(C) \ni r_j$ for some $g$ in
    the group, and $g(C)$ meets no earlier orbit (the orbits are invariant).
  - Prototype (`scratch/cc_bnb/`, `bb_blocks=2 bb_orbit=2`, the two orbits of the translations): double gross,
    $w = 15$: 9.76M nodes without symmetry, 1.15M with `start=0,144`, 0.58M with orbit exclusion; gross, $w = 12$:
    694K, 136K, 69K.
  - Not with `cbeg`/`cend`; codes only (see above for DEMs). `outC` holds the codewords found, not whole orbits (see
    the canonical representatives above); CC hits: one per codeword found.
  - Python: `orbits=` and `aut=` next to `start=`; the cache key `<key>:start=...` extended by a hash of the labels or
    of the generators; with verified generators, the result can be merged into the main record without `trust_start`.
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
- [ ] **Exact distance from extended CSS distance** Formulate explicit criteria
      on the exactness of the extended CSS distance, and learn to recognize such
      DEMs / circuits, as an alternative to code-sized model bounds
- [ ] **Optimize pair $p=2$ construction**: only start $i+j$ pairs with $i$ a verified codeword (non-trivial)
- [ ] **examine symmetry** of `model` check matrices to get guaranteed lower bounds faster (`--model-lb` option)
- [ ] **Branch-and-bound CC** (`cc_prune`), from G. Lee, "Bicycle Race: small two-block codes on the Pareto
  frontier," arXiv:2610.13880 (2026), App. C.4, whose exact engine uses the pruning of J. Y. Lee, K. Okada,
  N. Maskara, K. Kasai, H. Zhou, "Design principles for ultra-high-rate quantum codes," arXiv:2609.30069 (2026),
  App. B. Goal: `method=2` certification of distances up to $d = 20$ for matrices with row weights up to 12.
  - **Implemented in 0.15.0** (M2 of `scratch/cc_bnb/PLAN.md`): `cc_bnb()` in `dist_m4ri.c` as described below,
    `cc_prune=0|1|2` (default 1), the start-list exclusion, nodes per round with `debug&8` (exact: independent of the
    threads and work units; Python: `_last_run_stats["cc_rounds"]`), the Python pass-through (also
    `compute_model_lower_bound()`), Tests 110 and 111. M3 benchmark (2026-10-25): `benchmark/BENCHMARK.md` Section 9
    (4–35x for $d = 14$–18; the distances of the `yarn_*` codes with $d \ge 16$ and of `bb_360` certified, up to
    $d = 24$). Open: M4 (per-node cost, then the default of `cc_prune=2`); later runs: $d_X$ of `yarn_mitten_975`
    (the round $w = 23$, about 46 min), an independent check of the new distances with `cc_prune=0`, the start-list
    A/B, and the row-weight-12 GB code $[[254, 16, \le 24]]$.
  - First measurements (2026-10-25, 128 threads, idle server, one run each; `scratch/cc_bnb/m2/`; the node counts
    again with the exact counter, `nodes_count_w.sh`, as the nodes at the fixed levels of the work units were counted
    once per unit before): the same codewords as 0.14.8 in all 16 runs of `scratch/cc_bnb/m0/` (each codeword once);
    `method=2` up to the distance:
    - $X$ sectors: yarn_mitten_300 ($d = 14$) 1.69 s -> 0.23 s (round $w = 13$: 1.59G -> 49.8M nodes, 32x),
      yarn_ablp_560 ($d = 14$) 3.66 s -> 0.20 s (3.54G -> 47.5M nodes, 74x), double gross ($d = 18$) 1.36 s ->
      0.33 s (round $w = 17$: 15x fewer nodes), yarn_mitten_630 round $w = 11$: 151M -> 1.37M nodes (110x), 16.2 ->
      0.73 thread-s;
    - full DEMs: bb_90 ($d = 8$) 0.50 s -> 0.46 s; mitten_270 hef ($d = 7$) 5.5 s -> 2.6 s; mitten_150 custom
      ($d = 7$) 0.56 s -> 0.57 s (round $w = 6$: 1.5x fewer nodes, but 1.4x the CC work); `method=3`: gross circuit
      ($d = 10$) 189 s -> 79 s, mitten_270 custom ($d = 9$) 197 s -> 141 s;
    - `cc_prune=2`: 13–43% fewer nodes than 1, but 0.9–1.3x the time.
    - A node costs about 3x a node of `cc_walk()` (gross: 0.22 µs against 0.074 µs).
  - Prototype (Oct 2026, `scratch/cc_bnb/`: `bbcc.c`, one thread, with the objects of 0.14.7; plan in `PLAN.md`):
    `bb_F=0 bb_mrv=0 bb_b2=0` reproduces `cc_walk()` (gross, $w = 12$: 1884 codewords, 2084 hits). Lee's 133x and
    5000x were measured against an older dist-m4ri; the gains over 0.14.7 in the last round (node counts exact, times
    noisy on a loaded machine):
    - gross $w = 12$: 4x; double gross $w = 14$: 7–14x (45.4M -> 3.66M nodes); yarn_mitten_300 $w = 11$: 11–20x;
      yarn_smitten_330 $w = 11$: 17x; yarn_mitten_630_126_20, $w = 8, \dots, 11$: 8x, 15x, 28x, 64x (151M -> 1.1M
      nodes at $w = 11$); c204H: 2x.
    - Row weight 12 (`gen_gb.py`: generalized bicycle codes from weight-6 multiples of a primitive polynomial):
      $[[126,14,\le 14]]$, $w = 10$: 3.2x; $[[254,16,\le 24]]$, $w = 9$: 8x.
    - The growth per round drops from 5–7x to 3–4x, so that the gain grows with $w$; extrapolated to $w = 19$ for
      yarn_mitten_630: about $5\cdot 10^3$ s on one thread, against about $6\cdot 10^7$ s now.
    - Reduced DEMs (bb_72, bb_90, double gross): 0.86–1.2x for $w \le 6$, 1.6x at $w = 7$ (bb_90).
    - The relabeling of the columns and the order of the branches change the work by ±20%.
    - Contributions (single rounds, `scratch/cc_bnb/run_mrv.sh`, `mrv.log`; gross, double gross, yarn_mitten_300 and
      630, the GB code $[[254,16,\le 24]]$, the bb_72 and bb_90 DEMs): the fewest-available row alone (no $F$, no
      $b_2$) visits 1.05–1.7x fewer nodes, but takes 0.7–1.5x the time (as in `scratch/rowsort/`); adding $F$ gives
      1.1–1.5x fewer nodes again; $b_2$ alone (first row) 1.1–45x fewer nodes, 1.0–32x less time; all three 1.7–66x
      fewer nodes, 1.1–38x less time (the largest gains for the large $w$ of the mitten codes).
  - Kernel `cc_bnb()`, replacing `cc_walk()` (`smax=0`; `start_CC_recurs_mt()` stays for `smax>0`):
    - forbidden set $F$: branch $i$ adds the column $v_i$ and forbids $v_1, \dots, v_{i-1}$ in its subtree (also the
      columns whose branch was pruned or reached zero syndrome), so that every support is reached once;
    - branching row: the unsatisfied row with the fewest available columns (larger than the start column, not in
      $S \cup F$; `acnt[r]` updated in $O(\Delta)$ per column added or forbidden), a dead end if one has none;
      `sort_rows` then only breaks ties;
    - bounds with $R = w - |S|$ columns left: $b_1 = \lceil |s|/\Delta \rceil$ as now, and the packing bound $b_2$
      (unsatisfied rows with pairwise disjoint available sets, greedy with the smallest sets first, need distinct new
      columns): prune if $b_2 > R$; computed only for $|s| > R$ (otherwise it cannot prune) and $R \ge 3$ (the hash is
      cheaper for the last two levels);
    - the last two levels with the column-syndrome hash as now (with $F$: the second column never lies in the
      branching row), and the last column by a single lookup;
    - work units: at a fixed level, the branch `fix[l]` among the available columns of the row, with the earlier ones
      forbidden, so that the units still partition the tree;
    - `acnt` per worker, updated incrementally as the start column grows (a recount of $O(\mathrm{nnz})$ per start
      column made small rounds on the double-gross DEM 15x slower in the prototype).
  - `cc_prune=[0|1|2]`, default 1 for all inputs:
    - 0: the walk of 0.14.7 (first unsatisfied row in row order, $b_1$ only);
    - 1: $F$, the fewest-available row with the dead-end test, $b_1$ and $b_2$;
    - 2: also the coverage bound $b_3$ (the least $t$ such that the $t$ largest $c_v$ sum to at least $|s|$, with
      $c_v$ the number of unsatisfied rows of the available column $v$) and the fractional cover bound
      $\sum_r 1/\max_{v \in A_r} c_v$ over the unsatisfied rows $r$ (a dual solution of the covering LP): 30–60% fewer
      nodes than 1 at about the same time (±25%); the default if it wins after the per-node optimization;
    - Python: `cc_prune` passed like `sort_rows`, CLI `--cc-prune`.
  - `start` list: the search from a listed column forbids the columns listed before it (valid without further
    assumptions; the prototype shows a small gain, gross 136K -> 129K nodes, double gross 1.15M -> 1.11M; keep it if
    there is no slowdown), so that each codeword is found once.
  - CC hits: exactly one per codeword (needed for the RW miss estimate), as in all prototype runs with $F$; they enter
    `cnt` and `cw_hits_w` (the `min_hits` statistic) like `finC` codewords, not `cnt_u`/`cnt_w`. With `cc_prune=0`,
    codewords reached by several paths get more hits (gross: 2084 for 1884).
  - Per-node cost (a node costs about 2x a node of 0.14.7): one pass over the syndrome for the branching row and $b_2$,
    `acnt` as `uint16_t`, a tuned threshold for $b_2$, maybe $b_2$ updated from the packing of the parent.
  - `noscan=1` (the single round at `wmax`, certified with `dmin=wmax`) gets the same speed-up.
  - Tests: `cc_prune=0` and `1` give the same distances and the same irreducible codewords (`outC`) for the codes and
    DEMs of the benchmark, with `start`, `cbeg`/`cend`, `noscan`, `dmin`, `dW`, `maxC`, `method=3`, 1 or many threads;
    benchmark (1 and 128 threads): certify yarn_mitten_300, yarn_smitten_600_120_14, yarn_ablp_560, and try
    yarn_mitten_500, 540, 630.
  - Docs: README (CC section, `start`, help texts in `util_io.h`), `introduction.md`, Python docstrings; version 0.15.0
    in the four files.
- [x] **Irreducible codewords only** (codes and DEMs; 0.14.8, `irred=1` by default): a codeword which splits into a
  codeword and a zero-syndrome vector with disjoint supports is not recorded. The expert option `irred=0` of the binary
  records all codewords, as before (the Python CLI ignores it with a warning).
  - CC (`cc_record()`, both kernels): the columns of $H$ on the support $C$ have rank $|C| - 1$ (otherwise a proper
    subset of $C$ has zero syndrome); a small dense elimination on the at most $|C|\Delta$ rows touched
    (`csr_cols_rank()` in `util_m4ri.c`, with a per-thread workspace `cols_ws_t`), only for the codewords found.
  - RW pairs (`rw_pairs_step()`): the sum of two circuits is reducible iff their supports are disjoint, i.e., iff the
    weight of the sum is the sum of the weights (free); single RW columns give circuits, which are irreducible; the
    rows of the experimental `ksub` get the rank test.
  - `finC` (`nzlist_read()`, also in `dist_m4ri_old`): a reducible codeword is replaced by a nontrivial part (the
    lighter one if both parts are; `csr_cols_dependency()`, repeated until irreducible), which may lower `dmax`.
  - `debug&1`: the counts of the reducible codewords skipped by CC and RW, and of the `finC` codewords replaced.
  - Regression data (`scratch/cc_bnb/m0/`: `m0.sh`, the outputs of 0.14.7 for 16 deterministic runs of CC, RW, and
    `method=3` on codes and DEMs, and `verify_m1.sh`): the same results in all 16 runs; with `irred=0`, the same
    codeword sets as 0.14.7; with `irred=1`, exactly the reducible codewords removed (checked with `cwcheck`, which
    tests each codeword): surf_d3 DEM with `dW=3`: 102 of 1777 (CC), 848 of 2523 (RW with pairs, 3000 steps; none
    without pairs), 255 of 1930 (`method=3`), so that all three now give the same 1675 codewords; bb_72 DEM, RW with
    `dW=3`, 1000 steps: 15024 of 142284. None for the codes (CC at $w \le d + 2$) or at $w \le d + 2$ for the DEMs.
  - Before 0.14.8, there was no such filter: `nzlist_read()` checked the column range, the merged twins of a DEM,
    $Hc = 0$, and $Lc \neq 0$.
  - vecdec (`~/prog/vecdec/src/vecdec.c`) has two filters: `do_hash_decomp()` (after reading `finC`) splits each
    codeword into its connected components in the Tanner graph (`ufl_decompose()` in `dec_pre.c`), and replaces a
    reducible one (several components, or a trivial one) by its components with a nonzero syndrome of $L$;
    `do_hash_remove_reduc()` (after the RIS search) removes the codewords which contain a lighter stored codeword
    (at least `dx` = 1, or 3 for a DEM, lighter). CC codewords are connected by construction, so that the component
    test never fires on them: in the surf_d3 DEM, all the reducible codewords found by CC are connected (204, 2102,
    11552 hits at $w = 6, 7, 8$). The subset test finds them only if their lighter part is stored. Hence the rank
    test for CC and the RW pairs, and vecdec's replacement by the nontrivial parts for `finC`.
  - Minimum-weight codewords are irreducible; both kernels find every irreducible codeword (each proper prefix of its
    path has nonzero syndrome), but only some reducible ones, depending on the path: from $w = d + w_0$ on, with $w_0$
    the weight of the lightest zero-syndrome vector (3 for circuit DEMs).
  - Prototype (`bb_irr=1`): with the filter, both kernels give the same codewords at every weight (surf_d3 DEM,
    $w = 3, \dots, 8$), each hit once; without it, the sets differ from $w = 6$ on by reducible ones only (e.g., 204
    reducible hits at $w = 6$ for 0.14.7, 66 for the branch and bound).

- [ ] see if linear algebra can be used to eliminate the possibility of
      non-trivial DEM code words stretching over multiple REPEAT rounds.      
