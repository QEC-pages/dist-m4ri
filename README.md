# dist-m4ri - Distance of a Classical or Quantum CSS Code

## Overview

`dist-m4ri` is a high-performance multithreaded C program and Python library for computing and bracketing the minimum distance of binary classical linear codes, quantum CSS codes, and Stim Detector Error Models (DEMs).

The program implements three main methods:

- **Method 1 (`method=1`) - Random Window (RW) Algorithm**: Multithreaded random information set search to find low-weight non-trivial codewords and establish an **upper distance bound** $d_{\max}$.
- **Method 2 (`method=2`) - Connected Cluster (CC) Algorithm**: Multithreaded exhaustive depth-first cluster enumeration to compute **exact distance** or establish a certified **lower distance bound** $d_{\min}$.
- **Method 3 (`method=3`) - Bracketing Mode (Artillery Fork / Вилка)**: Concurrently runs CC and RW on multiple threads, dynamically balancing CPU cores between CC and RW based on current bounds $[d_{\min}, d_{\max}]$, distance estimate (`dexp`/`dest`), remaining RW steps, timeout, and the measured scaling characteristics of CC.

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
  - `dmin - 1` is the maximum cluster size analyzed without success by CC (`dmin = dmax` if CC found a minimum-weight codeword).
  - `dmax` is the weight of the smallest non-trivial codeword found by RW (`0` if none found).
  - `rw_steps` is the number of completed RW steps across all threads (`0` if CC found a minimum-weight codeword, or if RW did not run in `method=2`).
  - When `dmin = dmax = d`, the exact code distance is confirmed.

> **Note on Compatibility**: This 3-number output format (`dmin dmax rw_steps`) is specific to the multithreaded `dist_m4ri` and is incompatible with the legacy single-threaded `dist_m4ri_old` (which returned a single integer `d` or `-w`).

- **`stderr`**: Receives all status banners, thread balancing reports, CC round timings, RW discovery logs, warnings, and confinement profiles.

---

## How the Methods Work

### 1. Multithreaded RW Algorithm (`method=1`)
Searches for low-weight non-trivial binary codewords $c$ such that $Hc = 0$ and $Lc \neq 0$.
Threads independently generate random column permutations, compute Gaussian elimination to find information sets, and extract candidate dual-row codewords. When a lighter codeword is discovered, all worker threads atomically update the global upper bound $d_{\max}$ and prune heavier entries.

Relevant parameters:
- `steps=[int]`: Total number of information sets / RW rounds across all threads (default: 1).
- `wmin=[int]`: Minimum distance of interest (stop immediately when a codeword of weight $w \le w_{\min}$ is found).
- `threads=[int]`: Number of POSIX threads to run (default: number of CPU cores).
- `timeout=[sec]`: Maximum execution time in seconds (default: 60.0).

### 2. Multithreaded CC Algorithm (`method=2`)
Recursively explores connected clusters starting from each column $i \in [0, n-1]$. Columns are distributed dynamically among worker threads via lock-free atomic queues.
If `noscan=0` (default), CC scans weights $w = 1, 2, \dots, w_{\max}$. When `outC` is specified, CC exhausts all columns for weight $w$ to collect all unique minimum-weight codewords.

Relevant parameters:
- `wmax=[int]`: Maximum cluster weight to search (optional if `timeout>0` or `dmax>0` is specified; otherwise required for CC).
- `noscan=[int]`: If set to 1, start CC directly at $w_{\max}$ without scanning smaller weights.
- `cbeg=[int]`, `cend=[int]`: Column range $[c_{\text{beg}}, c_{\text{end}}]$ to limit the CC search space.
- `start=[int]`: Set $c_{\text{beg}} = c_{\text{end}} = \text{start}$ (useful for cyclic or symmetric codes).
- `smax=[int]`: Maximum syndrome weight to track for confinement profile (default: 5; set 0 to disable).

### 3. Bracketing Mode (`method=3`)
Dynamically partitions the available thread pool between CC (pushing $d_{\min}$ up) and RW (pulling $d_{\max}$ down) to solve the exact code distance as quickly as possible.

#### Dynamic Thread Allocation & Role of `dexp`:
1. **Target Search Depth**:
   - Before RW discovers a candidate codeword, $d_{\max}$ is unknown. Providing `dexp=D` (alias: `dest=D`) tells the coordinator to plan CC verification up to target weight $D - 1$ (or $D$).
2. **Predictive Workload Modeling**:
   - The coordinator measures empirical single-thread speed for RW steps ($t_{\text{RW}}$) and fits exponential growth to completed CC rounds to estimate time $T_{\text{CC}}$ required to reach $\min(d_{\max}, d_{\exp})$.
   - It computes the remaining work ratio:
     $$\text{ratio} = \frac{T_{\text{CC}}}{T_{\text{CC}} + T_{\text{RW}}}$$
     and dynamically splits threads at each round:
     $$N_{\text{CC}} = \text{round}(N_{\text{threads}} \times \text{ratio}), \quad N_{\text{RW}} = N_{\text{threads}} - N_{\text{CC}}$$
   - **Small $d_{\exp}$**: CC requires little work, so only 1–2 threads run CC while the majority maximize RW sampling speed.
   - **Heavier rounds**: As $w$ grows toward $d_{\exp}$, $T_{\text{CC}}$ increases and additional threads are shifted to CC to ensure both algorithms converge on the exact distance simultaneously.
3. **Adaptive Early Cutoff**:
   - If CC reaches $w > d_{\exp}$ before a codeword is found, CC halts and yields 100% of threads to RW.
   - If a projected CC round is estimated to exceed the remaining `timeout`, the coordinator terminates CC early and devotes remaining time entirely to RW.

Relevant parameters:
- `dexp=[int]` (alias: `dest=[int]`): Expected code distance to guide target search depth and thread allocation.
- `threads=[int]`: Number of worker threads (default: hardware concurrency).
- `timeout=[sec]`: Maximum execution time in seconds (default: 60.0).
- `steps=[int]`: Maximum total RW steps (default: 1000).
- `dW=[int]`: Extra weight window above $d_{\min}$ to continue collecting codewords ($w \le d_{\min} + \text{dW}$).

---

## Confinement Profile

With `smax > 0` (default: `smax=5`), the CC algorithm tracks the minimum non-zero syndrome weight observed for each cluster weight $w$:

```bash
$ ./src/dist_m4ri method=2 finH=./examples/surf_d5_H.mmx finL=./examples/surf_d5_L.mmx wmax=4 debug=0 threads=4
# confinement: 1,1,1,1
5 0
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

## Command-Line Usage

```sh
$ ./src/dist_m4ri --help
src/dist_m4ri: distance of a classical or quantum CSS code
	usage: src/dist_m4ri parameter=value [...]

   Required parameter:
	method=[int]: bitmap for method used (no default): 
		1: random window (RW) algorithm
		2: connected cluster (CC) algorithm
		3: bracketing mode (balanced concurrent RW and CC)

   Execution and multithreading parameters:
	threads=[int]: number of threads to use (0 for auto CPU count) (0)
	timeout=[sec]: timeout in seconds (60.0)
	dexp=[int]:    expected distance value for method=3 (alias: dest) (0)

   General parameters:
	finH=[str]: parity check matrix Hx (NULL)
	finG=[str]: matrix Hz (quantum CSS code only) (NULL)
	finL=[str]: matrix Lx (quantum CSS code only) (NULL)
	fdem=[str]: detector error model (DEM) file from stim (NULL)
	pmin=[float]: minimum error probability to keep for DEM (0.0)
	dmin=[int]: known lower bound on distance, inclusive (w starts from dmin in CC) (1)
	dmax=[int]: known upper bound on distance, inclusive (RW ignores codewords of weight >= dmax) (0)
	wmax=[int]: maximum cluster weight for CC (0; optional if timeout>0 or dmax>0)
	wmin=[int]: minimum distance of interest for RW and CC (1; early termination if cw <= wmin found)
	steps=[int]: number of RW information sets (1000)
	smax=[int]: maximum syndrome weight for confinement (5)
	outC=[str]: file to output codewords in NZLIST format (NULL)
	finC=[str]: file to read initial codewords from (NULL)
	dW=[int]: extra weight above dmin to collect codewords (-1)
	maxC=[int]: maximum number of codewords to collect (0)
	classical=[0|1]: force classical (1) or quantum (0) code (-1)
	debug=[int]: debug output bitmap (3)
```

### CLI Examples

```bash
# 1. Classical linear code using 8 threads in bracketing mode
$ ./src/dist_m4ri method=3 finH=./examples/c204H.mmx dest=10 steps=100000 threads=8
8 8

# 2. Stim Detector Error Model (DEM) with timeout and codeword export
$ ./src/dist_m4ri method=3 fdem=./examples/surf_d3.dem dexp=3 outC=cws.nz threads=4
3 3

# 3. Quantum CSS code (Hx and Hz) using pure CC search up to wmax=5
$ ./src/dist_m4ri method=2 finH=./examples/surf_d5_H.mmx finG=./examples/surf_d5_L.mmx wmax=5 threads=4
5 5
```

---

## Python Wrapper (`dist_m4ri.py`)

A Python module [`dist_m4ri.py`](dist_m4ri.py) is included for high-level scripting, NumPy/SciPy integration, and Stim interoperability without manual threading overhead.

### Key Python Functions

- `compute_classical_distance(H, ...)`: Minimum distance of a classical linear code (from NumPy 2D array, SciPy sparse matrix, or `.mtx` file).
- `compute_css_distance(Hx, Hz, Lx=None, Lz=None, ...)`: Distance $d = \min(d_X, d_Z)$ of a CSS quantum code.
- `compute_dem_distance(dem=None, circuit=None, ...)`: Minimum distance directly from a `stim.DetectorErrorModel`, `stim.Circuit`, or `.dem` file.
- `read_sparse_vectors(filepath)`: Parses NZLIST files into lists of 0-based integer support indices.
- Distance caching: `enable_distance_cache()`, `disable_distance_cache()`, `clear_distance_cache()`.
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
dist, d_info, cws = dist_m4ri.compute_dem_distance(circuit=circuit, do_cws=True, threads=8)
print(f"Surface code distance: {dist}, bounds: {d_info}, found {len(cws)} error mechanisms")

# 3. CSS Quantum Code (returns dist, dX_info, dZ_info where each info is [dmin, dmax, num_rw])
dist, d_x, d_z = dist_m4ri.compute_css_distance(
    Hx="examples/surf_d5_H.mmx",
    Hz="examples/surf_d5_H.mmx",
    Lz="examples/surf_d5_L.mmx",
    Lx="examples/surf_d5_L.mmx",
    d_exp=5,
    threads=8
)
print(f"CSS distance: {dist}, dX: {d_x}, dZ: {d_z}")  # CSS distance: 5, dX: [5, 5, 0], dZ: [5, 5, 0]

# 4. Persistent JSON Cache
# Pass cache_file to save/load bounds, cumulative RW steps, and codewords across runs:
d = dist_m4ri.compute_classical_distance("examples/c1920H.mmx", num_steps=100, cache_file="my_cache.json")
# Or set a default cache file globally:
dist_m4ri.set_distance_cache_file("my_cache.json")
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

# Run full C test suite (34 tests)
make test
```

### Python Unit Tests

```bash
pytest tests/test_dist_m4ri.py -v
```

---

## References

If you use this program, please cite:

*   A. Dumer, A. A. Kovalev, and L. P. Pryadko, "Distance verification for classical and quantum LDPC codes," *IEEE Transactions on Information Theory*, vol. 63, no. 7, pp. 4675-4690, 2017. [doi:10.1109/TIT.2017.2690381](https://doi.org/10.1109/TIT.2017.2690381).

Other related papers and software:

*   **vecdec Repository** (Random Information Set (RW) algorithm with error weights/probabilities):
    [QEC-pages/vecdec](https://github.com/QEC-pages/vecdec).

*   **QDistRnd GAP Package** (Random Information Set algorithm for quantum codes over arbitrary finite fields):
    L. P. Pryadko, V. A. Shabashov, and V. K. Kozin, "QDistRnd: A GAP package for computing the distance of quantum error-correcting codes," *Journal of Open Source Software*, vol. 7, no. 71, p. 4120, 2022. [doi:10.21105/joss.04120](https://doi.org/10.21105/joss.04120).

*   **Performance Comparison**:
    M. Webster, A. Jacob, and O. Higgott, "Distance-Finding Algorithms for Quantum Codes and Circuits," arXiv:2603.22532 [quant-ph], 2026. [arXiv:2603.22532](https://arxiv.org/abs/2603.22532).

