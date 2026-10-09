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

### Task 1: Subspace Sketching (`ksub`) for In-Cache Multicore Scaling
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
  1. Sample $k_{\text{sub}}$ rows uniformly from $N$ into the thread-local working matrix `M_sub`.
  2. Perform in-place RREF via `gauss_one(M_sub, perm->values[i], rank)` in permuted column order (stopping as soon as
     `rank == ksub`, avoiding matrix transpositions and non-thread-safe `m4ri_mmc` allocations).
  3. Directly extract candidate codewords from the reduced non-zero rows (already in sorted original column
     coordinates; no dual transposition, coordinate unpermutation, or `qsort` needed).
  4. Test logical non-triviality against $L$ (`sparse_syndrome_non_zero`), and atomically update `dmax`.

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
- [x] Add parameters `min_hits=[int]` (default `0` = disabled) and `cov_cws=[int]` (minimum distinct codewords
  required, default `1`).
- [x] **Codeword Multiplicity Tracking**:
  - Leverage the existing `cw_vec_t->cnt` field in `p->codewords` (which increments every time a duplicate codeword is
    found).
- [x] **Empirical Convergence Condition**:
  - In RW (`method=1` or `method=3`), track the number of distinct codewords found at the current minimum weight
    $d_{\text{min-found}}$ whose hit count reaches at least `min_hits`.
  - When at least `cov_cws` distinct minimum-weight codewords have each been independently rediscovered at least
    `min_hits` times (or all discovered minimum-weight codewords when `cov_cws <= 0`), signal early termination.
  - Optional variant: terminate when total repeat count or coverage ratio exceeds a statistical threshold, indicating
    that the minimum-weight codeword subspace has been saturated and the true code distance has likely been found.

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

## 2. Additional Algorithmic Explorations & Enhancements

- [ ] **Half-Weight Vector Meet-in-the-Middle**:
  - Investigate using the syndrome hash table or codeword hash to identify pairs of half-weight vectors whose sum forms
    a low-weight logical operator.
- [ ] **Trivial Codeword Statistics**:
  - Maintain counters of trivial codewords (stabilizers in quantum codes) of weights below $w_{\max}$ to monitor code
    degeneracy profile and ISD efficiency.
- [ ] **Sorting Optimization**:
  - Replace generic `qsort` on codeword support arrays `ee` with a small-array insertion sort or counting sort (since
    coordinates are already bounded by $n$ and array lengths are small, $w \le 100$).
- [ ] **Sparse Gaussian Elimination**:
  - Evaluate whether sparse elimination (e.g. CSR row combining) can outperform dense bit-matrices for very large,
    highly sparse DEMs where $k_{\text{sub}}$ is small.

---

## 3. Operational Modes Reference

### Mode: Distance Verification
- Given a suspected distance bound $d$ (e.g., $d = 10$), run with `method=1 wmax=10 wmin=9`.
- The program skips verification for codewords of weight $\ge 10$ and terminates immediately upon discovering any valid
  non-trivial codeword of weight $\le 9$ (outputting `-9`).
- If no codewords are found below $w_{\max}$, outputs `0`.

### Mode: Upper Bound Search with Hit-Count Saturation
- Run RW with `method=1 steps=N min_hits=M cov_cws=K`.
- Collects minimum-weight codewords into the hash table, terminating when the smallest observed weight has been
  confirmed by $M$ independent hits across $K$ distinct degenerate configurations.

### More debug bits
- Detailed timing information (measured / predicted values, reasons for termination)
- ???

