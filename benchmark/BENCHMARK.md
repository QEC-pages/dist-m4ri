# Benchmark Suite for `dist-m4ri`

This directory contains quantum LDPC CSS parity-check matrices (`.mtx`) and Stim circuit-level fault-tolerance files
(`.stim`) spanning code distances $d = 4$ to $d \le 34$. Detector Error Models (`.dem`)—both full `DEPOLARIZE2` DEMs
and simplified DEMs with minority-basis detectors stripped—can be generated on the fly from the `.stim` circuits using
`add_noise` in `dist_m4ri.py` (or `distance.py`).

> [!NOTE]
> The distances in the file names (and in the source patterns below) are the nominal values of the sources, and they
> are not always correct; the file names are kept unchanged. The distances computed with `dist_m4ri` are listed in the
> tables below: an exact value is certified (by CC), and $[d_{\min}, d_{\max}]$ is a bracket from a run which hit the
> timeout. The values of the Mitten codes and circuits are from a longer benchmark (`dist_m4ri` 0.12.1, 128 threads,
> up to 2.8 h per run); the other values are from `method=3` runs on 16 threads (30–60 s per run), and they are
> provisional.

## Directory Structure

- `css/`: MatrixMarket (`.mtx`) CSS parity-check matrices (`_Hx.mtx`, `_Hz.mtx`) and canonical logical operators
  (`_Lx.mtx`, `_Lz.mtx` for Mitten codes).
- `stim_bb/`: Unmodified Bivariate Bicycle (BB), Gross, and Double-Gross `.stim` circuits from external repositories.
- `stim_mitten/`: Mitten code `.stim` circuits generated via `quits` with `CircuitBuildOptions(get_all_detectors=True)`.
- `stim_torus/`: 2D local torus Bivariate Bicycle `.stim` circuits (`d = 14` to `17`, noiseless templates).
- `figures/`: scaling of `dist_m4ri` with the number of threads on some of these codes and circuits (Section 4).

---

## Origins and References

1. **Bivariate Bicycle (BB), Gross (`[[144, 12, 12]]`), and Double-Gross (`[[288, 12, 18]]`) Codes**
   - **Code Construction**: S. Bravyi, A. W. Cross, J. M. Gambetta, D. Maslov, P. Rall, and T. J. Yoder,
     *High-threshold and low-overhead fault-tolerant quantum memory*, Nature **627**, 778–782 (2024), arXiv:2308.07915.
   - **Higher-Distance BB Codes (`[[360, 12, 24]]`, `[[756, 16, <=34]]`) & `si1000` Stim Circuits**:
     Copied unmodified from [`quantumlib/tesseract-decoder`](https://github.com/quantumlib/tesseract-decoder)
     (commit `e7c762eef24161e304ba6fccb856a05b41f88c39`, directory `testdata/bivariatebicyclecodes/`), associated
     with A. Beni et al., arXiv:2503.10988. These circuits use the superconducting-inspired `si1000` noise model at
     $p = 0.0001$ with $r = d$ rounds in both $X$ and $Z$ memory bases.
   - **Uniform-Noise (`uniform_circuit`) BB Stim Circuits**:
     Copied unmodified from [`trmue/relay`](https://github.com/trmue/relay)
     (commit `d185194ba0cb4101ced4340d82b2ee6d42f225f0`, directory `tests/testdata/bicycle_bivariate/`), associated
     with T. Müller et al., arXiv:2503.11514. These circuits use standard uniform circuit-level depolarizing noise at
     $p = 0.001$ with $r = d$ rounds in the $X$ memory basis.

2. **Mitten Codes (Non-Abelian Lifted-Product CSS Codes)**
   - **Code Construction & Stim Circuits**: Generated via [`mkangquantum/quits`](https://github.com/mkangquantum/quits)
     (commit `0f36b8c5a9b81d6b7be0cef2c81c61b5850c0452`), associated with *Quantum LDPC codes for modular
     architectures via lifted products over non-abelian groups*, arXiv:2607.28795.
   - Constructed as $1 \times 2$ lifted products $LP(A, B)$ over non-abelian groups $G = C_m \times S_3$ of order
     $|G| = 6m$, yielding parameters $[[5|G|, |G|, d]]$ with row weight 9 and column weight 3 or 6.
   - Circuits are generated at $p = 0.001$ (`ErrorModel(1e-3, 1e-3, 1e-3, 1e-3)`, `basis="Z"`,
     `CircuitBuildOptions(get_all_detectors=True)`) for two syndrome-extraction schedules:
     - `custom`: depth-12 interleaved mixed $X/Z$ schedule.
     - `hef` (`hook_error_free`): depth-24 group-element-layer schedule where each entangling layer is a permutation
       matching between one check block and one data block.

3. **2D Local Torus Bivariate Bicycle Stim Circuits**
   - Copied unmodified from `tmp/*.stim`. These are noiseless syndrome-extraction templates on a 2D torus with Swap/CNOT
     routing (`L3_torus`) for codes with the nominal parameters $[[60, 8, 14]]$, $[[120, 8, 15]]$, $[[200, 16, 16]]$,
     and $[[200, 8, 17]]$ (from the file names; the code distances have not been verified) in both $X$ and $Z$ bases.
     When passed to `dist_m4ri.compute_dem_distance(circuit=...)` or
     `dist_m4ri.add_noise(circuit, p=0.001)`, standard circuit-level depolarizing noise (`DEPOLARIZE1`, `DEPOLARIZE2`,
     `X_ERROR`/`Z_ERROR`) is added on the fly.

---

## 1. CSS Parity-Check and Logical Matrices (`benchmark/css/`)

| Prefix | $[[n, k, d]]$ | Size ($m \times n$) | Construction / Polynomials or Group |
| :--- | :--- | :--- | :--- |
| `bb_72_12_6` | $[[72, 12, 6]]$ | $36 \times 72$ | $\ell=6, m=6$, $A=x^3+y+y^2$, $B=y^3+x+x^2$ |
| `bb_90_8_10` | $[[90, 8, 10]]$ | $45 \times 90$ | $\ell=15, m=3$, $A=x^9+y+y^2$, $B=1+x^2+x^7$ |
| `bb_108_8_10` | $[[108, 8, 10]]$ | $54 \times 108$ | $\ell=9, m=6$, $A=x^3+y+y^2$, $B=y^3+x+x^2$ |
| `gross_144_12_12` | $[[144, 12, 12]]$ | $72 \times 144$ | $\ell=12, m=6$, $A=x^3+y+y^2$, $B=y^3+x+x^2$ |
| `double_gross_288_12_18` | $[[288, 12, 18]]$ | $144 \times 288$ | $\ell=12, m=12$, $A=x^3+y^2+y^7$, $B=y^3+x+x^2$ |
| `bb_360_12_24` | $[[360, 12, \le 24]]$ | $180 \times 360$ | $\ell=30, m=6$, $A=x^9+y+y^2$, $B=y^3+x^{25}+x^{26}$ |
| `bb_756_16_34` | $[[756, 16, \le 34]]$ | $378 \times 756$ | $\ell=21, m=18, A=x^3+y^{10}+y^{17}, B=y^5+x^3+x^{19}$ |
| `mitten_90_18_6_C3xS3` | $[[90, 18, 6]]$ | $36 \times 90$ | $G = C_3 \times S_3$ ($|G|=18$, `seed=2`) |
| `mitten_150_30_8_C5xS3` | $[[150, 30, 8]]$ | $60 \times 150$ | $G = C_5 \times S_3$ ($|G|=30$, `seed=2`) |
| `mitten_210_42_8_C7xS3` | $[[210, 42, 6]]$ | $84 \times 210$ | $G = C_7 \times S_3$ ($|G|=42$, `seed=2`) |
| `mitten_270_54_10_C9xS3` | $[[270, 54, 10]]$ | $108 \times 270$ | $G = C_9 \times S_3$ ($|G|=54$, `seed=2`) |
| `mitten_330_66_6_C11xS3` | $[[330, 66, 4]]$ | $132 \times 330$ | $G = C_{11}\times S_3$ ($|G|=66$, `seed=2`) |

Distances computed with `dist_m4ri`:

| Prefix | $d_X$ | $d_Z$ |
| :--- | :--- | :--- |
| `bb_72_12_6` | 6 | 6 |
| `bb_90_8_10` | 10 | 10 |
| `bb_108_8_10` | – | 10 |
| `gross_144_12_12` | 12 | 12 |
| `double_gross_288_12_18` | 18 | 18 |
| `bb_360_12_24` | – | [18, 24] |
| `bb_756_16_34` | – | [17, 34] |
| `mitten_90_18_6_C3xS3` | 6 | 6 |
| `mitten_150_30_8_C5xS3` | 8 | 8 |
| `mitten_210_42_8_C7xS3` | **6** | 8 |
| `mitten_270_54_10_C9xS3` | 10 | 10 |
| `mitten_330_66_6_C11xS3` | **4** | 6 |

- $d_X$ ($d_Z$) is the minimum weight of an $X$-type ($Z$-type) logical operator, i.e., of a vector $c$ with
  $H_Z c = 0$ ($H_X c = 0$) which is not in the row space of $H_X$ ($H_Z$); $d = \min(d_X, d_Z)$. "–": not computed
  yet.
- For the BB codes, $d_X = d_Z$ by the symmetry exchanging the two sectors, so one sector suffices.
- `mitten_330_66_6_C11xS3` has distance **4**, not 6: $d_X = 4$ (132 $X$-type logical operators of weight 4 were
  found, and those checked independently over GF(2) are valid), while $d_Z = 6$.
- `mitten_210_42_8_C7xS3` has distance **6**, not 8: $d_X = 6$, while $d_Z = 8$.
- The Mitten distances are certified by CC (under 0.5 s per sector on 128 threads), and the lightest logical operator
  found in each run passed an independent check over GF(2). Numbers of minimum-weight logical operators ($X$ / $Z$,
  enumerated by CC): 66 / 90 (`mitten_90`), 180 / 120 (`mitten_150`), 84 / 42 (`mitten_210`), 216 / 648
  (`mitten_270`), and 132 / 132 (`mitten_330`).
- For `bb_360` and `bb_756`, the published distances are upper bounds (the codewords found by `dist_m4ri` have the
  same weights), and the lower bounds are from 60 s runs.

---

## 2. Stim Circuits (`stim_bb/`, `stim_mitten/`, `stim_torus/`)

The circuit-level distances are listed in [Section 2.4](#24-circuit-level-distances).

### 2.1 Bivariate Bicycle, Gross, and Double-Gross (`benchmark/stim_bb/`)

| File (`benchmark/stim_bb/`) | Repo | Noise | $p$ | Rounds | Basis | Source Pattern |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `bb_72_12_6_si1000_r6_X.stim` | `tesseract` | `si1000` | `1e-4` | 6 | `X` | `r=6,d=6,nkd=[[72,12,6]]` |
| `bb_72_12_6_si1000_r6_Z.stim` | `tesseract` | `si1000` | `1e-4` | 6 | `Z` | `r=6,d=6,nkd=[[72,12,6]]` |
| `bb_90_8_10_si1000_r10_X.stim` | `tesseract` | `si1000` | `1e-4` | 10 | `X` | `r=10,d=10,nkd=[[90,8,10]]` |
| `bb_90_8_10_si1000_r10_Z.stim` | `tesseract` | `si1000` | `1e-4` | 10 | `Z` | `r=10,d=10,nkd=[[90,8,10]]` |
| `bb_108_8_10_si1000_r10_X.stim` | `tesseract` | `si1000` | `1e-4` | 10 | `X` | `r=10,d=10,nkd=[[108,8,10]]` |
| `bb_108_8_10_si1000_r10_Z.stim` | `tesseract` | `si1000` | `1e-4` | 10 | `Z` | `r=10,d=10,nkd=[[108,8,10]]` |
| `gross_144_12_12_si1000_r12_X.stim` | `tesseract` | `si1000` | `1e-4` | 12 | `X` | `r=12,d=12,nkd=[[144,12,12]]` |
| `gross_144_12_12_si1000_r12_Z.stim` | `tesseract` | `si1000` | `1e-4` | 12 | `Z` | `r=12,d=12,nkd=[[144,12,12]]` |
| `double_gross_288_12_18_si1000_r18_X.stim` | `tesseract` | `si1000` | `1e-4` | 18 | `X` | `r=18,d=18,[[288,12,18]]` |
| `double_gross_288_12_18_si1000_r18_Z.stim` | `tesseract` | `si1000` | `1e-4` | 18 | `Z` | `r=18,d=18,[[288,12,18]]` |
| `bb_360_12_24_si1000_r24_X.stim` | `tesseract` | `si1000` | `1e-4` | 24 | `X` | `r=24,d=24,nkd=[[360,12,24]]` |
| `bb_360_12_24_si1000_r24_Z.stim` | `tesseract` | `si1000` | `1e-4` | 24 | `Z` | `r=24,d=24,nkd=[[360,12,24]]` |
| `bb_756_16_34_si1000_r34_X.stim` | `tesseract` | `si1000` | `1e-4` | 34 | `X` | `r=34,d=34,nkd=[[756,16,34]]` |
| `bb_756_16_34_si1000_r34_Z.stim` | `tesseract` | `si1000` | `1e-4` | 34 | `Z` | `r=34,d=34,nkd=[[756,16,34]]` |
| `bb_72_12_6_uniform_r6_X.stim` | `relay` | `uniform` | `1e-3` | 6 | `X` | `..._72_12_6_memory_X` |
| `gross_144_12_12_uniform_r12_X.stim` | `relay` | `uniform` | `1e-3` | 12 | `X` | `..._144_12_12_memory_X` |
| `double_gross_288_12_18_uniform_r18_X.stim` | `relay` | `uniform` | `1e-3` | 18 | `X` | `..._288_12_18_memory_X` |

### 2.2 Mitten Codes (`benchmark/stim_mitten/`)

| File (`benchmark/stim_mitten/`) | Source | Schedule | Depth | $p$ | Rounds | Basis |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `mitten_90_18_6_C3xS3_custom_r4_Z.stim` | `quits` | `custom` (interleaved) | 12 | `1e-3` | 4 | `Z` |
| `mitten_90_18_6_C3xS3_hef_r4_Z.stim` | `quits` | `hook_error_free` | 24 | `1e-3` | 4 | `Z` |
| `mitten_150_30_8_C5xS3_custom_r5_Z.stim` | `quits` | `custom` (interleaved) | 12 | `1e-3` | 5 | `Z` |
| `mitten_150_30_8_C5xS3_hef_r5_Z.stim` | `quits` | `hook_error_free` | 24 | `1e-3` | 5 | `Z` |
| `mitten_210_42_8_C7xS3_custom_r4_Z.stim` | `quits` | `custom` (interleaved) | 12 | `1e-3` | 4 | `Z` |
| `mitten_210_42_8_C7xS3_hef_r4_Z.stim` | `quits` | `hook_error_free` | 24 | `1e-3` | 4 | `Z` |
| `mitten_270_54_10_C9xS3_custom_r4_Z.stim` | `quits` | `custom` (interleaved) | 12 | `1e-3` | 4 | `Z` |
| `mitten_270_54_10_C9xS3_hef_r4_Z.stim` | `quits` | `hook_error_free` | 24 | `1e-3` | 4 | `Z` |
| `mitten_330_66_6_C11xS3_custom_r4_Z.stim` | `quits` | `custom` (interleaved) | 12 | `1e-3` | 4 | `Z` |
| `mitten_330_66_6_C11xS3_hef_r4_Z.stim` | `quits` | `hook_error_free` | 24 | `1e-3` | 4 | `Z` |

The Mitten circuits are $Z$-basis memory experiments, whose logical observables are flipped by $X$-type errors; their
circuit-level distance cannot exceed the code's $d_X$, which is 6 for `mitten_210` and 4 for `mitten_330` (Section 1).

### 2.3 2D Local Torus BB Circuits (`benchmark/stim_torus/`)

| File (`benchmark/stim_torus/`) | Nominal $[[n, k, d]]$ | Basis | Notes |
| :--- | :--- | :--- | :--- |
| `NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_X.stim` | $[[60, 8, 14]]$ | `X` | Noiseless (`tmp/`) |
| `NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_Z.stim` | $[[60, 8, 14]]$ | `Z` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w12p_n120_k8_d15_X.stim` | $[[120, 8, 15]]$ | `X` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w12p_n120_k8_d15_Z.stim` | $[[120, 8, 15]]$ | `Z` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k16_d16_X.stim` | $[[200, 16, 16]]$ | `X` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k16_d16_Z.stim` | $[[200, 16, 16]]$ | `Z` | Noiseless (`tmp/`) |
| `ENENESESENENE_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k8_d17_X.stim` | $[[200, 8, 17]]$ | `X` | Noiseless (`tmp/`) |
| `ENENESESENENE_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k8_d17_Z.stim` | $[[200, 8, 17]]$ | `Z` | Noiseless (`tmp/`) |

The nominal code parameters are taken from the file names; the code distances of the torus codes have not been
verified yet (no check matrices are included).

### 2.4 Circuit-Level Distances

The circuit-level distance is the minimum number of faults (error mechanisms of the DEM) which cause an undetectable
logical error.
- It is computed for the stripped DEM (minority-basis detectors removed), which is the default of
  `dist_m4ri.compute_dem_distance(circuit=...)`, with 30 s `method=3` runs on 16 threads (Mitten circuits: up to 2.1 h
  on 128 threads).
- In parentheses: the distance of the full DEM, where it was computed. The stripped DEM keeps a subset of the
  detectors of the full DEM, so its distance is a lower bound on the full-DEM distance. †: lower bound from the
  stripped DEM, upper bound from a codeword of the full DEM (Mitten circuits: found by RW on the full DEM in 240 s on
  128 threads; CC on the full DEMs of `mitten_270`, for 2.8 h, did not close the brackets).
- "–": not computed yet; "n/a": no circuit for this basis.
- For the noiseless torus templates, the noise of `add_noise(p=0.001)` is added. The circuit distance depends on the
  noise model (which faults are included), but not on $p$.

| Circuit (`.stim` file name without `_X`/`_Z`) | $X$ basis | $Z$ basis |
| :--- | :--- | :--- |
| `bb_72_12_6_si1000_r6` | 6 (full: 6†) | 6 |
| `bb_90_8_10_si1000_r10` | [7, 8] | – |
| `bb_108_8_10_si1000_r10` | [7, 8] | – |
| `gross_144_12_12_si1000_r12` | [7, 10] | – |
| `double_gross_288_12_18_si1000_r18` | [4, 24] (full: [4, 53]†) | – |
| `bb_360_12_24_si1000_r24` | [4, 29] | – |
| `bb_756_16_34_si1000_r34` | – | – |
| `bb_72_12_6_uniform_r6` | 6 | n/a |
| `gross_144_12_12_uniform_r12` | [7, 11] | n/a |
| `double_gross_288_12_18_uniform_r18` | – | n/a |
| `mitten_90_18_6_C3xS3_custom_r4` | n/a | 5 (full: 5) |
| `mitten_90_18_6_C3xS3_hef_r4` | n/a | 4 (full: 4) |
| `mitten_150_30_8_C5xS3_custom_r5` | n/a | 7 (full: 7†) |
| `mitten_150_30_8_C5xS3_hef_r5` | n/a | 5 (full: 5†) |
| `mitten_210_42_8_C7xS3_custom_r4` | n/a | 6 (full: 6†) |
| `mitten_210_42_8_C7xS3_hef_r4` | n/a | 5 (full: 5†) |
| `mitten_270_54_10_C9xS3_custom_r4` | n/a | 9 (full: [9, 10]†) |
| `mitten_270_54_10_C9xS3_hef_r4` | n/a | 7 (full: [7, 8]†) |
| `mitten_330_66_6_C11xS3_custom_r4` | n/a | 4 (full: 4†) |
| `mitten_330_66_6_C11xS3_hef_r4` | n/a | 4 (full: 4†) |
| `NNEESEEESEENN_L3_torus_..._n60_k8_d14` | 6 (full: 6†) | 6 |
| `NEENNWNWNNEEN_L3_torus_..._n120_k8_d15` | 6 | – |
| `NEENNWNWNNEEN_L3_torus_..._n200_k16_d16` | [6, 8] | – |
| `ENENESESENENE_L3_torus_..._n200_k8_d17` | [6, 8] | – |

The circuit-level distances are often below the code distances in the file names. For example, the gross-code
circuit has distance at most 10, consistent with the circuit-level distance $\le 10$ reported for the IBM
syndrome-extraction circuit. For the Mitten codes, the "hook-error-free" schedule has a lower circuit-level distance
than the `custom` schedule (4 vs 5, 5 vs 7, 5 vs 6, and 7 vs 9 for `mitten_90`, `mitten_150`, `mitten_210`, and
`mitten_270`; 4 with both schedules for `mitten_330`). The `custom` schedule reaches the code's $d_X$ for `mitten_210`
(6) and `mitten_330` (4), and $d_X - 1$ for `mitten_90`, `mitten_150`, and `mitten_270` (unless the full DEM of
`mitten_270` has distance 10).

---

## 3. Generating Full vs. Minority-Detector-Stripped DEMs on the Fly

In circuit-level noise models using `DEPOLARIZE2(p)`, cross-term two-qubit Pauli errors (such as $X \otimes Z$,
$Y \otimes X$, $Y \otimes Z$, $Y \otimes Y$) occur with probability $p/15 = O(p)$, the same order as single-type errors
($X \otimes X$ or $Z \otimes Z$). Consequently, the full DEM (`circuit.detector_error_model(decompose_errors=False)`)
couples $X$-detectors and $Z$-detectors into a single large parity-check matrix with up to $\sim 7\times$ more unique
error columns than a single-basis projection.

To avoid storing 41 MB of pre-generated `.dem` files in the repository, DEMs can be generated on the fly from the
`.stim` files using `dist_m4ri.add_noise` (for noiseless circuits) and `strip_minority_detectors` (in `distance.py`):
```python
import stim
import dist_m4ri

circ = stim.Circuit.from_file("benchmark/stim_torus/NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_X.stim")
if not dist_m4ri.has_noise(circ):
    circ = dist_m4ri.add_noise(circ, p=0.001)
dem_full = circ.detector_error_model(decompose_errors=False)
dist, bounds = dist_m4ri.compute_dem_distance(dem=dem_full)
```

| Stim Circuit (`.stim` -> `.dem`) | Basis | Full (`dets x errs`) | Stripped (`dets x errs`) | $k$ |
| :--- | :--- | :--- | :--- | :--- |
| `mitten_90_18_6_C3xS3_custom_r4_Z.dem` | `Z` | $360 \times 23,394$ | $216 \times 3,330$ | 18 |
| `mitten_90_18_6_C3xS3_hef_r4_Z.dem` | `Z` | $360 \times 22,410$ | $216 \times 3,330$ | 18 |
| `mitten_150_30_8_C5xS3_custom_r5_Z.dem` | `Z` | $720 \times 47,700$ | $420 \times 6,660$ | 30 |
| `mitten_150_30_8_C5xS3_hef_r5_Z.dem` | `Z` | $720 \times 46,130$ | $420 \times 6,660$ | 30 |
| `mitten_210_42_8_C7xS3_custom_r4_Z.dem` | `Z` | $840 \times 54,726$ | $504 \times 7,770$ | 42 |
| `mitten_210_42_8_C7xS3_hef_r4_Z.dem` | `Z` | $840 \times 52,437$ | $504 \times 7,665$ | 42 |
| `mitten_270_54_10_C9xS3_custom_r4_Z.dem` | `Z` | $1,080 \times 70,254$ | $648 \times 9,990$ | 54 |
| `mitten_270_54_10_C9xS3_hef_r4_Z.dem` | `Z` | $1,080 \times 67,482$ | $648 \times 9,990$ | 54 |
| `mitten_330_66_6_C11xS3_custom_r4_Z.dem` | `Z` | $1,320 \times 85,932$ | $792 \times 12,210$ | 66 |
| `mitten_330_66_6_C11xS3_hef_r4_Z.dem` | `Z` | $1,320 \times 81,642$ | $792 \times 11,880$ | 66 |
| `bb_72_12_6_si1000_r6_X.dem` | `X` | $432 \times 17,280$ | $252 \times 2,664$ | 12 |
| `bb_72_12_6_uniform_r6_X.dem` | `X` | $432 \times 16,200$ | $252 \times 2,232$ | 12 |
| `bb_90_8_10_si1000_r10_X.dem` | `X` | $900 \times 37,440$ | $495 \times 5,490$ | 8 |
| `bb_108_8_10_si1000_r10_X.dem` | `X` | $1,080 \times 44,928$ | $594 \times 6,588$ | 8 |
| `gross_144_12_12_si1000_r12_X.dem` | `X` | $1,728 \times 72,576$ | $936 \times 10,512$ | 12 |
| `gross_144_12_12_uniform_r12_X.dem` | `X` | $1,728 \times 67,824$ | $936 \times 8,784$ | 12 |
| `double_gross_288_12_18_si1000_r18_X.dem` | `X` | $5,184 \times 221,184$ | $2,736 \times 31,392$ | 12 |
| `double_gross_288_12_18_uniform_r18_X.dem` | `X` | $5,184 \times 206,496$ | $2,736 \times 26,208$ | 12 |
| `bb_360_12_24_si1000_r24_X.dem` | `X` | $8,640 \times 371,520$ | $4,500 \times 52,200$ | 12 |
| `bb_756_16_34_si1000_r34_X.dem` | `X` | $25,704 \times 1,112,832$ | $13,230 \times 154,980$ | 16 |
| `NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_X.dem` | `X` | $120 \times 4,830$ | $90 \times 1,500$ | 8 |
| `NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_Z.dem` | `Z` | $120 \times 4,830$ | $90 \times 1,500$ | 8 |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w12p_n120_k8_d15_X.dem` | `X` | $240 \times 9,600$ | $180 \times 3,000$ | 8 |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w12p_n120_k8_d15_Z.dem` | `Z` | $240 \times 9,600$ | $180 \times 3,000$ | 8 |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k16_d16_X.dem` | `X` | $400\times 16,000$ | $300\times 5,000$ |16|
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k16_d16_Z.dem` | `Z` | $400\times 16,000$ | $300\times 5,000$ |16|
| `ENENESESENENE_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k8_d17_X.dem` | `X` | $400 \times 14,800$ | $300 \times 4,800$ | 8|
| `ENENESESENENE_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k8_d17_Z.dem` | `Z` | $400 \times 14,800$ | $300 \times 4,800$ | 8|

---

## 4. Scaling with the Number of Threads (`benchmark/figures/`)

Wall time and speedup $T(1)/T(p)$ of `dist_m4ri` for $p = 1, \dots, 128$ threads on a two-socket AMD EPYC 7B13 server
(2 x 32 cores, 128 logical CPUs with SMT), with `nothrottle=1`, `threads=p`, and `taskset` to the logical CPUs
$0..p-1$ (socket 0 first, then socket 1, then the SMT siblings for $p > 64$; the gray line marks $p = 64$, one thread
per physical core). Each run does a fixed amount of work: `method=1` a fixed number of RW `steps` with `min_hits=0`,
`method=2` all CC rounds up to `wmax`, and `method=3` both (or, with the default parameters, the search for the exact
distance). Versions: 0.12.2 (`method=1` only for $p = 16, 64, 128$: its RW kernel is that of 0.12.1), 0.12.1, the
master branch (0.9.0), and the legacy single-threaded `dist_m4ri_old` (one thread). The times include reading the
input (median of repeated runs for the fast cases).

![Scaling for the gross code](figures/scaling_gross_css.png)

**gross** $[[144, 12, 12]]$, one sector (`css/gross_144_12_12_Hx.mtx` with `_Hz.mtx`; `method=1`: `steps=3000000`,
`method=2`: `wmax=11`, `method=3`: default parameters). RW scales almost linearly up to 64 threads (58x), and the SMT
siblings add another 40% (82x at 128 threads). CC up to $w = 11$ takes only 0.22 s on one thread: beyond 8–16
threads, the run (25–45 ms) is limited by the start-up of the threads. Compared with master, RW is 1.9x faster, and CC
21x (one thread).

![Scaling for the double-gross code](figures/scaling_dgross_css.png)

**double-gross** $[[288, 12, 18]]$, one sector (`steps=1000000`; `wmax=15`; `steps=500000 wmax=15`). RW scales to 90x
at 128 threads. In 0.12.1, CC stops at a speedup of 16: a CC round hands out start columns, and the clusters from the
first start column, which can grow to all larger columns, are about $w/n \approx 5\%$ of the round. With the CC work
units of 0.12.2 (start columns split into pieces), CC up to $w = 15$ takes 0.44 s on 128 threads (speedup 89) instead
of 3.1 s, and `method=3` 0.75 s instead of 3.2 s; master needs 485 s for this CC run on 128 threads.

![Scaling for the double-gross circuit DEM](figures/scaling_dgross_dem.png)

**double-gross circuit** (`stim_bb/double_gross_288_12_18_si1000_r18_X.stim`, stripped DEM $2736 \times 31392$;
`steps=300`; `wmax=5`; `steps=150 wmax=5`). RW on this large matrix is limited by the memory bandwidth (dense copies
of $H$ and $H^T$, about 21 MB per thread): at most about 20x, at 64 threads on both sockets. Pinned to one socket, 16
and 32 threads take the same time, and the SMT siblings make RW slower; with the default parameters (no pinning), the
memory throttle limits RW to 32 threads on this DEM. CC scales well (48x at 128 threads) and is 20–30x faster than
master.
