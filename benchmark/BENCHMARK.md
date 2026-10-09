# Benchmark Suite for `dist-m4ri`

This directory contains quantum LDPC CSS parity-check matrices (`.mtx`) and Stim circuit-level fault-tolerance files
(`.stim`) spanning code distances $d = 4$ to $d \le 34$. Detector Error Models (`.dem`)—both full `DEPOLARIZE2` DEMs
and simplified DEMs with minority-basis detectors stripped—can be generated on the fly from the `.stim` circuits using
`add_noise` in `dist_m4ri.py` (or `distance.py`).

> [!NOTE]
> The distances in the file names (and in the source patterns below) are the nominal values of the sources, and they
> are not always correct; the file names are kept unchanged. The distances computed with `dist_m4ri` (version 0.10.1,
> `method=3`, 16 threads, 30–60 s per run) are listed in the tables below: an exact value is certified (by CC), and
> $[d_{\min}, d_{\max}]$ is a bracket from a run which hit the timeout. These numbers are provisional and will be
> updated after a longer benchmark.

## Directory Structure

- `css/`: MatrixMarket (`.mtx`) CSS parity-check matrices (`_Hx.mtx`, `_Hz.mtx`) and canonical logical operators
  (`_Lx.mtx`, `_Lz.mtx` for Mitten codes).
- `stim_bb/`: Unmodified Bivariate Bicycle (BB), Gross, and Double-Gross `.stim` circuits from external repositories.
- `stim_mitten/`: Mitten code `.stim` circuits generated via `quits` with `CircuitBuildOptions(get_all_detectors=True)`.
- `stim_torus/`: 2D local torus Bivariate Bicycle `.stim` circuits (`d = 14` to `17`, noiseless templates).

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
     routing (`L3_torus`) for $[[60, 8, 14]]$, $[[120, 8, 15]]$, $[[200, 16, 16]]$, and $[[200, 8, 17]]$ codes in both
     $X$ and $Z$ bases. When passed to `dist_m4ri.compute_dem_distance(circuit=...)` or
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
| `bb_360_12_24` | $[[360, 12, 24]]$ | $180 \times 360$ | $\ell=30, m=6$, $A=x^9+y+y^2$, $B=y^3+x^{25}+x^{26}$ |
| `bb_756_16_34` | $[[756, 16, \le 34]]$ | $378 \times 756$ | $\ell=21, m=18, A=x^3+y^{10}+y^{17}, B=y^5+x^3+x^{19}$ |
| `mitten_90_18_6_C3xS3` | $[[90, 18, 6]]$ | $36 \times 90$ | $G = C_3 \times S_3$ ($|G|=18$, `seed=2`) |
| `mitten_150_30_8_C5xS3` | $[[150, 30, 8]]$ | $60 \times 150$ | $G = C_5 \times S_3$ ($|G|=30$, `seed=2`) |
| `mitten_210_42_8_C7xS3` | $[[210, 42, 8]]$ | $84 \times 210$ | $G = C_7 \times S_3$ ($|G|=42$, `seed=2`) |
| `mitten_270_54_10_C9xS3` | $[[270, 54, 10]]$ | $108 \times 270$ | $G = C_9 \times S_3$ ($|G|=54$, `seed=2`) |
| `mitten_330_66_6_C11xS3` | $[[330, 66, 6]]$ | $132 \times 330$ | $G = C_{11}\times S_3$ ($|G|=66$, `seed=2`) |

---

## 2. Stim Circuits (`stim_bb/`, `stim_mitten/`, `stim_torus/`)

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

### 2.3 2D Local Torus BB Circuits (`benchmark/stim_torus/`)

| File (`benchmark/stim_torus/`) | $[[n, k, d]]$ | Basis | Notes |
| :--- | :--- | :--- | :--- |
| `NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_X.stim` | $[[60, 8, 14]]$ | `X` | Noiseless (`tmp/`) |
| `NNEESEEESEENN_L3_torus_s1_NW_w10p_s2_SW_w6p_n60_k8_d14_Z.stim` | $[[60, 8, 14]]$ | `Z` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w12p_n120_k8_d15_X.stim` | $[[120, 8, 15]]$ | `X` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w12p_n120_k8_d15_Z.stim` | $[[120, 8, 15]]$ | `Z` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k16_d16_X.stim` | $[[200, 16, 16]]$ | `X` | Noiseless (`tmp/`) |
| `NEENNWNWNNEEN_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k16_d16_Z.stim` | $[[200, 16, 16]]$ | `Z` | Noiseless (`tmp/`) |
| `ENENESESENENE_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k8_d17_X.stim` | $[[200, 8, 17]]$ | `X` | Noiseless (`tmp/`) |
| `ENENESESENENE_L3_torus_s1_NW_w10p_s2_SW_w20p_n200_k8_d17_Z.stim` | $[[200, 8, 17]]$ | `Z` | Noiseless (`tmp/`) |

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
