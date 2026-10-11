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
  - Separate hit statistics for the uniform and the windowed RW steps.
  - Empirical per-weight hit rates from the representative set, extrapolated to weights below the minimum found.
  - Account for the small-weight trivial vectors explicitly (see **Trivial Codeword Statistics** below).
  - To check: the `ksub` basis refresh (Task 3) puts found low-weight codewords into the basis $N$, which may inflate
    their hit counts and $\langle n\rangle$.

---

## 2. Additional Algorithmic Explorations & Enhancements

- [ ] **Half-Weight Vector Meet-in-the-Middle**:
  - Investigate using the syndrome hash table or codeword hash to identify pairs of half-weight vectors whose sum forms
    a low-weight logical operator.
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
- [ ] **Distance Benchmark Update**:
  - Run the longer benchmark to certify or tighten the provisional distances in `benchmark/BENCHMARK.md` (plan in
    `tmp/benchmark_distance_plan.md`).

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
  - [ ] Python could parse the `# code parameters:` line instead of computing ranks itself.
- [x] in addition, print estimated prob to find (miss) a single codeword of given
  weight based on the number of information sets and the n and distance found with RW.
  - Done: Python `--verbose`, and the binary with the hit uniformity warning (Task 4); a proper estimate for DEMs
    is still open (Task 6).
- [ ] The periodic status shows that in `method=3` a CC round can wait for a whole RW batch before any CC thread
  starts (e.g., `c1920H.mmx`: 100 RW steps of 6.6 ms per batch, i.e., 0.66 s): smaller RW batches while CC rounds
  are short, or an RW batch interruption when a CC round starts.

### More improvements 
- [ ] try verifying (using cluster decomposition from vecdec (?)) if the
      constructed / imported non-trivial CW can be decomposed into disconnected
      components.  Keep only non-trivial part (does it ever happen? how
      expensive is this?)
- [ ] See that matrices constructed when generating L from H (e.g.) are reused when codeword searching.
- [ ] For a very large circuit, when RW is very slow, should we try constructing
      codewords from the underlying code (can this even be done from a DEM ---
      or should it be done in python?)
- [ ] 
