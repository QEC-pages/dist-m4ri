"""
Unit tests for dist_m4ri.py Python wrapper.
"""

import json
import os
import re
import shutil
import sys
from pathlib import Path
import pytest
import numpy as np
import scipy.sparse as sp

# Add root directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import dist_m4ri

EXAMPLES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "examples"))


def test_find_binary():
    bin_path = dist_m4ri.find_dist_m4ri_binary()
    assert os.path.isfile(bin_path)
    assert os.access(bin_path, os.X_OK)


def test_classical_distance_file():
    h_file = os.path.join(EXAMPLES_DIR, "c204H.mmx")
    d = dist_m4ri.compute_classical_distance(h_file, d_exp=10, threads=4)
    assert d == 8


def test_classical_distance_numpy():
    # Hamming [7, 4, 3] code
    H = np.array([
        [1, 0, 0, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 1]
    ], dtype=np.int8)
    d = dist_m4ri.compute_classical_distance(H, threads=2)
    assert d == 3

    # With codewords
    d, cws = dist_m4ri.compute_classical_distance(H, do_cws=True, threads=2)
    assert d == 3
    assert len(cws) > 0
    assert all(len(cw) == 3 for cw in cws)


def test_css_distance_files():
    hx_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    hz_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    # For surf_d5: Hx as finH, Hz as finL gives d=5
    dist, d_x, d_z = dist_m4ri.compute_css_distance(
        Hx=hx_file, Hz=hx_file, Lz=hz_file, Lx=hz_file,
        d_exp=5, threads=4
    )
    assert dist == 5
    assert d_x == [5, 5, 0]
    assert d_z == [5, 5, 0]


def test_css_distance_sparse():
    # Surface code d=3
    hx = sp.csr_matrix([
        [1, 1, 0, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 1, 1, 0],
        [0, 0, 0, 0, 1, 1, 0, 1, 1]
    ], dtype=np.int8)
    hz = sp.csr_matrix([
        [1, 0, 0, 1, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 1, 0, 0, 1]
    ], dtype=np.int8)
    dist, d_x, d_z = dist_m4ri.compute_css_distance(hx, hz, threads=2)
    assert dist > 0
    assert len(d_x) == 3
    assert len(d_z) == 3


def test_dem_distance_file():
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    dist, d_info = dist_m4ri.compute_dem_distance(dem=dem_file, threads=4)
    assert dist == 3
    assert d_info == [3, 3, 0]

    # With codewords (method=3 bracketing mode collects discovered cws)
    dist, d_info, cws = dist_m4ri.compute_dem_distance(dem=dem_file, do_cws=True, threads=4)
    assert dist == 3
    assert d_info == [3, 3, 0]
    assert len(cws) > 0
    assert all(len(cw) == 3 for cw in cws)

    # Exhaustive CC scan (method=2) finds all 128 codewords with drop=0 (a fresh cache: the codewords of this run)
    dist_m4ri.clear_distance_cache()
    dist, d_info, cws_cc = dist_m4ri.compute_dem_distance(dem=dem_file, method=2, wmax=3, do_cws=True, threads=4,
                                                          drop=0)
    assert dist == 3
    assert d_info == [3, 3, 0]
    assert len(cws_cc) == 128
    assert all(len(cw) == 3 for cw in cws_cc)

    # With drop=auto (the default), the binary drops 8 of the 24 detectors: the 21 codewords of the reduced DEM (its
    # 55 columns touch no dropped detector) are codewords of the DEM
    dist_m4ri.clear_distance_cache()
    dist, d_info, cws_red = dist_m4ri.compute_dem_distance(dem=dem_file, method=2, wmax=3, do_cws=True, threads=4)
    assert (dist, d_info, len(cws_red)) == (3, [3, 3, 0], 21)
    assert set(map(tuple, cws_red)) <= set(map(tuple, cws_cc))


def test_dem_distance_stim():
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=3,
        distance=3,
        after_clifford_depolarization=0.001
    )
    dem = circuit.detector_error_model(decompose_errors=True)
    dist, d_info = dist_m4ri.compute_dem_distance(dem=dem, threads=4)
    assert dist == 3
    assert d_info == [3, 3, 0]


def test_caching():
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()

    H = np.array([
        [1, 0, 0, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 1]
    ], dtype=np.int8)

    d1 = dist_m4ri.compute_classical_distance(H, threads=2)
    assert len(dist_m4ri._distance_cache) == 1

    # Second call should be a cache hit
    d2 = dist_m4ri.compute_classical_distance(H, threads=2)
    assert d1 == d2
    assert len(dist_m4ri._distance_cache) == 1

    dist_m4ri.clear_distance_cache()
    assert len(dist_m4ri._distance_cache) == 0


def test_codedistance_solver():
    if not dist_m4ri._HAS_CODEDISTANCE:
        pytest.skip("codedistance is not installed")

    H = np.array([
        [1, 0, 0, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 1]
    ], dtype=np.int8)

def test_run_dist_m4ri_three_numbers():
    # Method 2 (CC only): rw_steps must be 0
    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    dmin, dmax, rw_steps = dist_m4ri.run_dist_m4ri(method=2, finH=h_file, finL=l_file, wmax=5, threads=4)
    assert (dmin, dmax, rw_steps) == (5, 5, 0)

    # Method 2 (CC not found): dmin=wmax+1, dmax=0, rw_steps=0
    dmin, dmax, rw_steps = dist_m4ri.run_dist_m4ri(method=2, finH=h_file, finL=l_file, wmax=3, threads=4)
    assert (dmin, dmax, rw_steps) == (4, 0, 0)

    # Method 1 (RW): rw_steps reported (min_hits=0: all steps; the reduced DEM, drop=auto, has only 55 columns)
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    dmin, dmax, rw_steps = dist_m4ri.run_dist_m4ri(method=1, fdem=dem_file, steps=100, threads=4, min_hits=0)
    assert dmin == 1
    assert dmax == 3
    assert rw_steps >= 100


def test_dmin_dmax_parameters():
    # Test dmin/dmax in run_dist_m4ri
    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    dmin, dmax, rw_steps = dist_m4ri.run_dist_m4ri(
        method=3, finH=h_file, finL=l_file, dmin=4, dmax=5, timeout=5, threads=4
    )
    assert (dmin, dmax) == (5, 5)

    # Test dmin/dmax in compute_classical_distance
    c_file = os.path.join(EXAMPLES_DIR, "c204H.mmx")
    d = dist_m4ri.compute_classical_distance(c_file, dmin=5, dmax=8, threads=4)
    assert d == 8

    # Test dmin/dmax in compute_dem_distance
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    d_dem, _ = dist_m4ri.compute_dem_distance(dem=dem_file, dmin=2, dmax=4, threads=4)
    assert d_dem == 3


def test_caching_cumulative_rw_steps():
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()

    c_file = os.path.join(EXAMPLES_DIR, "c1920H.mmx")

    # Run 1: 50 steps
    d1 = dist_m4ri.compute_classical_distance(c_file, method=1, num_steps=50, threads=4)
    entry1 = dist_m4ri.get_cached_distance(H=c_file)
    assert entry1 is not None
    assert entry1["rw_steps"] == 50
    assert entry1["dmax"] > 0
    assert d1 == entry1["dmax"]

    # Run 2: another 100 steps
    d2 = dist_m4ri.compute_classical_distance(c_file, method=1, num_steps=100, threads=4)
    entry2 = dist_m4ri.get_cached_distance(H=c_file)
    assert entry2 is not None
    assert entry2["rw_steps"] == 150
    assert entry2["dmax"] <= entry1["dmax"]
    assert d2 == entry2["dmax"]

    dist_m4ri.clear_distance_cache()


def test_format_bounds_list_and_str():
    # Exact known distance
    assert dist_m4ri.format_bounds_list(5, 5, 0) == [5, 5, 0]
    assert dist_m4ri.format_bounds_str([5, 5, 0]) == "5 5 0 (exact)"

    # No upper bound (dmax == 0)
    assert dist_m4ri.format_bounds_list(4, 0, 0) == [4, 0, 0]
    assert dist_m4ri.format_bounds_str([4, 0, 0]) == "4 0 0"

    # No lower bound (dmin <= 1)
    assert dist_m4ri.format_bounds_list(0, 323, 100) == [0, 323, 100]
    assert dist_m4ri.format_bounds_str([0, 323, 100]) == "0 323 100"

    # Lower and upper bounds differing
    assert dist_m4ri.format_bounds_list(4, 6, 1000) == [4, 6, 1000]
    assert dist_m4ri.format_bounds_str([4, 6, 1000]) == "4 6 1000"


def test_persistent_json_cache(tmp_path):
    import json
    json_file = str(tmp_path / "test_cache.json")

    dist_m4ri.clear_distance_cache()
    c_file = os.path.join(EXAMPLES_DIR, "c1920H.mmx")

    # Run 1: 50 steps with cache_file
    d1 = dist_m4ri.compute_classical_distance(c_file, method=1, num_steps=50, threads=4, cache_file=json_file)
    assert os.path.isfile(json_file)

    with open(json_file, "r") as f:
        data1 = json.load(f)
    assert data1.get("__version__") == dist_m4ri.__version__
    cache_entries1 = {k: v for k, v in data1.items() if not k.startswith("__")}
    assert len(cache_entries1) == 1
    key = list(cache_entries1.keys())[0]
    assert data1[key]["rw_steps"] == 50
    assert data1[key]["dmax"] == d1

    # Clear memory cache and re-read from JSON
    dist_m4ri.clear_distance_cache()
    entry = dist_m4ri.get_cached_distance(H=c_file, cache_file=json_file)
    assert entry is not None
    assert entry["rw_steps"] == 50

    # Run 2: another 100 steps
    d2 = dist_m4ri.compute_classical_distance(c_file, method=1, num_steps=100, threads=4, cache_file=json_file)
    with open(json_file, "r") as f:
        data2 = json.load(f)
    assert data2[key]["rw_steps"] == 150
    assert data2[key]["dmax"] <= d1

    # CSS persistent caching
    hx_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    hz_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    # First run with method=2 up to wmax=3 (lower bound dmin=4, dmax=0): must NOT be treated as exact
    d_css_lb, dx_lb, dz_lb = dist_m4ri.compute_css_distance(
        Hx=hx_file, Hz=hx_file, Lz=hz_file, Lx=hz_file,
        method=2, wmax=3, threads=4, cache_file=json_file
    )
    assert d_css_lb == 4
    assert dx_lb == [4, 0, 0]
    assert dz_lb == [4, 0, 0]
    entry_css_lb = dist_m4ri.get_cached_distance(Hx=hx_file, Hz=hx_file, Lx=hz_file, Lz=hz_file, cache_file=json_file)
    assert entry_css_lb["dmin"] == 4
    assert entry_css_lb["dmax"] == 0

    # Second run with method=1 (RW): must execute RW (not falsely hit exact cache) and merge bounds [4, 5, steps]
    d_css_rw, dx_rw, dz_rw = dist_m4ri.compute_css_distance(
        Hx=hx_file, Hz=hx_file, Lz=hz_file, Lx=hz_file,
        method=1, num_steps=50, min_hits=0, threads=4, cache_file=json_file
    )
    assert d_css_rw == 5
    assert dx_rw[0] == 4 and dx_rw[1] == 5 and dx_rw[2] >= 50
    assert dz_rw[0] == 4 and dz_rw[1] == 5 and dz_rw[2] >= 50

    # Third run with method=3: certifies exact d=5 -> [5, 5, 0]
    d_css, dx_info, dz_info = dist_m4ri.compute_css_distance(
        Hx=hx_file, Hz=hx_file, Lz=hz_file, Lx=hz_file,
        d_exp=5, threads=4, cache_file=json_file
    )
    assert d_css == 5
    assert dx_info == [5, 5, 0]
    assert dz_info == [5, 5, 0]
    with open(json_file, "r") as f:
        data_css = json.load(f)
    assert len(data_css) >= 2

    # DEM persistent caching
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    d_dem, d_info = dist_m4ri.compute_dem_distance(dem=dem_file, threads=4, cache_file=json_file)
    assert d_dem == 3
    with open(json_file, "r") as f:
        data_dem = json.load(f)
    assert len(data_dem) >= 3

    dist_m4ri.clear_distance_cache(cache_file=json_file, clear_file=True)
    assert not os.path.exists(json_file)


def test_quantum_distance_single_sided():
    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    dist, d_info = dist_m4ri.compute_quantum_distance(
        H=h_file, L=l_file, method=3, d_exp=5, threads=4, return_info=True
    )
    assert dist == 5
    assert d_info == [5, 5, 0]


def test_cli_argument_parsing():
    # Auto-infer classical = 0 when finG or finL is given
    args1 = dist_m4ri.parse_cli_args(["finH=h.mtx", "finG=g.mtx", "smax=0", "finC=init.nz", "start=2"])
    assert args1["classical"] == 0
    assert args1["finH"] == "h.mtx"
    assert args1["finG"] == "g.mtx"
    assert args1["smax"] == 0
    assert args1["finC"] == "init.nz"
    assert args1["start"] == [2]  # start is a list of CC start columns (since version 0.10.0)

    # Auto-infer classical = 1 when only finH is given
    args2 = dist_m4ri.parse_cli_args(["finH=h.mtx", "method=2"])
    assert args2["classical"] == 1
    assert args2["finH"] == "h.mtx"
    assert args2["finG"] is None

    # Auto-infer classical = 0 when fdem is given
    args3 = dist_m4ri.parse_cli_args(["fdem=model.dem"])
    assert args3["classical"] == 0


def test_debug_bits(capsys):
    # debug values: decimal or hexadecimal (no octal), OR-combined in the CLI
    assert dist_m4ri.parse_debug_value("0x10004") == dist_m4ri.PY_DBG_COMMANDS | 4
    assert dist_m4ri.parse_debug_value("011") == 11
    for bad in ("-1", "0x", "abc", "1.5"):
        with pytest.raises(ValueError):
            dist_m4ri.parse_debug_value(bad)
    args = dist_m4ri.parse_cli_args(["fdem=model.dem", "debug=1", "debug=0x10000"])
    assert args["debug"] == 1 | dist_m4ri.PY_DBG_COMMANDS
    # the Python bits are not passed to the binary: the command line is echoed with debug=0 (binary silent)
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    res = dist_m4ri.run_dist_m4ri(method=2, fdem=dem_file, wmax=3, threads=2, debug=dist_m4ri.PY_DBG_COMMANDS)
    assert res == (3, 3, 0)
    out = capsys.readouterr().out
    assert " debug=0 method=2 " in out and "Finished in" in out
    assert "# input:" not in out
    # binary summary bit: the input summary (the DEM reduced by drop=auto, the default) and why the run ended
    dist_m4ri.run_dist_m4ri(method=2, fdem=dem_file, wmax=3, threads=2, debug=1)
    out = capsys.readouterr().out
    assert "# input: detector error model, n=55" in out
    assert "dropped 8 of 24 detectors (drop=auto), exact" in out
    assert "CC found min-weight codeword: d=3" in out
    dist_m4ri.run_dist_m4ri(method=2, fdem=dem_file, wmax=3, threads=2, debug=1, drop="0")
    assert "# input: detector error model, n=221" in capsys.readouterr().out
    # keep (and list) the temporary files
    dist_m4ri.clear_distance_cache()
    dist_m4ri.disable_distance_cache()
    try:
        H = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.int8)
        d = dist_m4ri.compute_classical_distance(H, method=2, wmax=3, threads=1, debug=dist_m4ri.PY_DBG_KEEP_FILES)
        assert d == 3
        out = capsys.readouterr().out
        kept = [line.split("Kept temporary file: ", 1)[1] for line in out.splitlines() if "Kept temporary" in line]
        assert kept and all(os.path.exists(f) for f in kept)
        for f in kept:
            os.remove(f)
    finally:
        dist_m4ri.enable_distance_cache()


def test_quantum_cache_separation():
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()

    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")

    # Classical distance of H
    d_class = dist_m4ri.compute_classical_distance(h_file, dmax=5, threads=4)
    # Quantum distance of (H, L)
    d_quant = dist_m4ri.compute_quantum_distance(h_file, L=l_file, dmax=5, threads=4)

    cache = dist_m4ri.get_distance_cache()
    class_keys = [k for k in cache if k.startswith("classical:")]
    quant_keys = [k for k in cache if k.startswith("quantum:")]

    assert len(class_keys) >= 1
    assert len(quant_keys) >= 1

    dist_m4ri.clear_distance_cache()


def test_explain_bounds():
    # Exact distance
    exp1 = dist_m4ri.explain_bounds([5, 5, 0], method=2, label="dX")
    assert "Lower bound (dmin = 5): Exact distance certified" in exp1
    assert "Upper bound (dmax = 5): Weight of the smallest non-trivial codeword discovered" in exp1
    assert "Random window steps (rw_steps = 0): Set to 0 because the exact distance d = 5 was proven" in exp1

    # Pure lower bound in method 2
    exp2 = dist_m4ri.explain_bounds([4, 0, 0], method=2)
    assert "All cluster weights w <= 3 were exhaustively analyzed" in exp2
    assert "Method 2 (Connected Cluster) is an exhaustive search" in exp2

    # RW search
    exp3 = dist_m4ri.explain_bounds([1, 8, 120], method=1)
    assert "No non-trivial lower bound certified" in exp3
    assert "120 completed random information set searches" in exp3


def test_rw_infoset_estimate():
    # P1(w) = w C(n-w, k-1) / C(n, k) with k = n - rank: surf_d3.dem (n=221, rank(H)=24), w=3
    from math import comb
    p1 = dist_m4ri.rw_infoset_find_prob(221, 24, 3)
    assert p1 == pytest.approx(3 * comb(218, 196) / comb(221, 197), rel=1e-9)
    assert dist_m4ri.rw_infoset_find_prob(221, 24, 26) == 0.0  # w > rank + 1: never found
    # c204 (n=204, rank(H)=101, d=8): the hardest lighter weight is w=7; miss probability (1-P1)^steps
    est = dist_m4ri.rw_infoset_estimate(204, 101, 6000, 8, step_time=1e-5, steps_total=12000)
    p7 = dist_m4ri.rw_infoset_find_prob(204, 101, 7)
    assert est["w"] == 7
    assert est["p_find"] == pytest.approx(p7)
    assert est["p_miss"] == pytest.approx((1.0 - p7) ** 6000)
    n1 = int(np.ceil(np.log(0.01) / np.log1p(-p7)))
    assert est["steps_target"] == int(np.ceil(n1 * 2.0))  # half of the steps uniform
    assert est["time_target"] == pytest.approx(est["steps_target"] * 1e-5)
    # weights below a certified dmin are excluded; nothing to check if dmin >= d
    assert dist_m4ri.rw_infoset_estimate(204, 101, 6000, 8, dmin=8) is None
    # explanation from the stderr statistics of the binary
    stderr = (
        "# codewords accumulated: total=100, min_w=8: cws=1, total_hits=170, hits min=170, max=170, avg=170.00, "
        "stdev=0.00, <n>=5.27 over 100 cws of w=8..15 (min_hits=5, cov_cws=100)\n"
        "# RW information sets: n=204, rank(H)=101, steps=6000 (uniform permutations: 6000), ksub=0, "
        "4.4e-05 s/step per thread, 8 RW threads\n"
    )
    stats = dist_m4ri._parse_stderr_stats(stderr)
    assert stats["rw_n"] == 204 and stats["rw_rank"] == 101 and stats["rw_steps_uniform"] == 6000
    assert stats["rw_pairs_l"] == 0  # (the summary line of a binary without the RW pairs)
    exp = dist_m4ri.explain_bounds([0, 8, 6000], method=1, stats=stats)
    assert "Hit-based estimate: P_fail ~ exp(-<n>) = 0.0051 (below 1%)" in exp
    assert "a single codeword of weight 7 < 8 is found with probability" in exp
    assert "Model check: expected hits per codeword of weight 8" in exp
    # with the RW pairs (l pivot rows of the buckets): P = P1 + P2, P2 = C(w,2) C(n-w,k-2) / C(n,k) C(r-w+2,l) / C(r,l)
    p2 = comb(7, 2) * comb(197, 101) / comb(204, 103) * comb(96, 4) / comb(101, 4)
    assert dist_m4ri.rw_infoset_find_prob(204, 101, 7, pairs_l=4) == pytest.approx(p7 + p2, rel=1e-9)
    assert dist_m4ri.rw_infoset_find_prob(204, 101, 7, pairs_l=0) == pytest.approx(p7, rel=1e-12)
    est2 = dist_m4ri.rw_infoset_estimate(204, 101, 6000, 8, pairs_l=4)
    assert est2["w"] == 7 and est2["p_find"] == pytest.approx(p7 + p2)
    stderr2 = stderr.replace("ksub=0, ", "ksub=0, pairs=4, ") + (
        "# RW pairs (p=2, buckets of l=4 random pivot rows): 1.23e+03 pairs checked per step; codewords recorded with "
        "weight up to the upper bound: 517 from pairs, 170 from single columns\n")
    stats2 = dist_m4ri._parse_stderr_stats(stderr2)
    assert stats2["rw_pairs_l"] == 4 and stats2["rw_step_time"] == pytest.approx(4.4e-05)
    assert stats2["rw_threads"] == 8
    assert (stats2["rw_pairs_per_step"], stats2["rw_hits_pairs"], stats2["rw_hits_single"]) == (1230.0, 517, 170)
    exp2 = dist_m4ri.explain_bounds([0, 8, 6000], method=1, stats=stats2)
    assert "Information-set estimate (uniform random information sets with pairs, l = 4; n = 204" in exp2
    assert f"is found with probability {p7 + p2:.3g} per RW step" in exp2
    assert ("RW pairs (p=2, buckets of l = 4 pivot rows): 1.23e+03 pairs checked per RW step; codewords recorded "
            "with weight up to the upper bound: 517 from pairs, 170 from single columns.") in exp2


def test_cli_css_dx_dz_bounds(capsys):
    hx_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    lz_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    ret = dist_m4ri.main([
        f"Hx={hx_file}", f"Hz={hx_file}", f"Lx={lz_file}", f"Lz={lz_file}",
        "method=2", "wmax=5", "--no-cache", "threads=4"
    ])
    assert ret == 0
    captured = capsys.readouterr()
    assert "dX: 5 5 0 (exact)" in captured.out
    assert "dZ: 5 5 0 (exact)" in captured.out
    assert "(d = 5) (exact)" in captured.out


def test_cli_verbose_mode(capsys):
    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    ret = dist_m4ri.main([
        "--verbose", "--no-cache", "method=3", f"finH={h_file}", f"finL={l_file}",
        "wmax=5", "threads=4", "min_hits=3", "cov_cws=5"
    ])
    assert ret == 0
    captured = capsys.readouterr()
    assert "Cache retrieval:" in captured.out
    assert "Lower bound" in captured.out
    assert "Upper bound" in captured.out
    assert "Random window steps" in captured.out
    assert "Codewords accumulated:" in captured.out
    assert "hits min =" in captured.out
    assert "avg =" in captured.out
    assert "stdev =" in captured.out


def test_check_finc_outc():
    import tempfile

    # Non-existent file when finC == outC
    assert dist_m4ri.check_finc_outc("nonexistent.nz", "nonexistent.nz") is None

    # Non-existent file when finC != outC
    assert dist_m4ri.check_finc_outc("nonexistent.nz", "other.nz") == "nonexistent.nz"

    # Empty 0-byte file when finC == outC
    with tempfile.NamedTemporaryFile(suffix=".nz", delete=False) as f:
        tmp_empty = f.name
    try:
        assert dist_m4ri.check_finc_outc(tmp_empty, tmp_empty) is None
    finally:
        if os.path.exists(tmp_empty):
            os.remove(tmp_empty)

    # Existing non-empty file when finC == outC
    with tempfile.NamedTemporaryFile(suffix=".nz", delete=False, mode="w") as f:
        f.write("1 2 3\n")
        tmp_nonempty = f.name
    try:
        assert dist_m4ri.check_finc_outc(tmp_nonempty, tmp_nonempty) == tmp_nonempty
    finally:
        if os.path.exists(tmp_nonempty):
            os.remove(tmp_nonempty)


def test_cli_identical_finc_outc_nonexistent(capsys):
    import tempfile
    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    
    tmp_out = os.path.join(tempfile.gettempdir(), f"tmp_test_cw_{os.getpid()}.nz")
    if os.path.exists(tmp_out):
        os.remove(tmp_out)

    try:
        ret = dist_m4ri.main([
            "--verbose", "--no-cache", "method=2",
            f"finH={h_file}", f"finL={l_file}",
            f"finC={tmp_out}", f"outC={tmp_out}",
            "wmax=5", "threads=4"
        ])
        assert ret == 0
        captured = capsys.readouterr()
        assert "Warning: finC=" in captured.out
        assert "is empty or non-existent; silently ignoring input codewords." in captured.out
        assert "5 5 0 (exact)" in captured.out
        assert os.path.exists(tmp_out)
    finally:
        if os.path.exists(tmp_out):
            os.remove(tmp_out)


def test_cli_help(capsys):
    ret = dist_m4ri.main(["--help"])
    assert ret == 0
    captured = capsys.readouterr()
    assert "--morehelp" in captured.out
    assert "threads=N" in captured.out
    # Check that it fits within 80 rows
    lines = [line for line in captured.out.splitlines() if line.strip()]
    assert len(lines) < 80


def test_cli_morehelp(capsys):
    ret = dist_m4ri.main(["--morehelp"])
    assert ret == 0
    captured = capsys.readouterr()
    assert "Required input" in captured.out
    assert "Calculation method:" in captured.out
    assert "threads=N" in captured.out
    assert "nothrottle=1" in captured.out


def test_cli_morehelp_topic(capsys):
    def entries(out):  # (section headers, the first words of the entries)
        lines = out.splitlines()
        return ([s for s in lines if s and not s[0].isspace()],
                [s.split()[0] for s in lines if s.startswith("  ") and not s.startswith("   ")])

    # --morehelp TOPIC: the entries of TOPIC, under their section header, with the other entries of its group (see
    # also) and the group (section 'Option groups')
    groups = "Option groups (--help GROUP prints all entries of GROUP):"
    assert dist_m4ri.main(["--morehelp", "maxC"]) == 0
    assert entries(capsys.readouterr().out) == (
        ["Codeword collection and export:", groups],
        ["--cws", "outC=FILE", "finC=FILE", "maxC=N", "dW=N", "irred=0|1", "expand=0|1", "codewords"])
    # several topics (case-insensitive, aliases): in the order of the full help, each section header once
    assert dist_m4ri.main(["--morehelp=SEED,method,v"]) == 0
    assert entries(capsys.readouterr().out) == (["Calculation method:", "General options:"],
                                                ["method=1|2|3", "--verbose", "seed=N"])
    # '-', '_', and '=VALUE' are ignored; aliases in the descriptions ('Alias: batch=N', 'also trust_start=1')
    for topic, first in (("no-cache", "--no-cache"), ("--rounds=3", "--rounds"), ("batch", "chunk_size=N"),
                         ("trust_start", "--trust-start"), ("v", "--verbose"), ("out_dem=1", "--out-dem")):
        assert dist_m4ri.main([f"--morehelp={topic}"]) == 0
        assert first in entries(capsys.readouterr().out)[1], topic
    # a group (also singular or plural): all its entries; a parameter of a group gives the whole group
    for topic in ("windows", "window", "kwin", "WIN_HYB"):
        assert dist_m4ri.main(["-h", topic]) == 0
        assert entries(capsys.readouterr().out) == (["Multithreading & throttling:", groups],
                                                    ["kwin=N", "win_mode=0|1", "win_hyb=N", "windows"]), topic
    assert dist_m4ri.main(["--help", "experimental"]) == 0
    assert entries(capsys.readouterr().out)[1] == ["ksub=N", "refresh=N", "experimental"]
    assert dist_m4ri.main(["--help", "windos"]) == 1
    assert "no help entry for 'windos'; did you mean windows" in capsys.readouterr().err
    # an unknown topic: the closest names on stderr, exit code 1 (the known topics are still printed)
    assert dist_m4ri.main(["--morehelp", "mthod"]) == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "no help entry for 'mthod'; did you mean method?" in captured.err
    assert dist_m4ri.main(["--morehelp", "maxC,xyz"]) == 1
    captured = capsys.readouterr()
    assert "maxC=N" in entries(captured.out)[1]
    assert "no help entry for 'xyz'" in captured.err
    # an option after --morehelp is not a topic: the full help
    assert dist_m4ri.main(["--morehelp", "--verbose"]) == 0
    assert "Calculation method:" in capsys.readouterr().out
    # --help TOPIC (also -h TOPIC, --help=TOPIC, help=TOPIC): the same as --morehelp TOPIC; topics are collected
    assert dist_m4ri.main(["--morehelp", "maxC,dest"]) == 0
    expected = capsys.readouterr().out
    assert {"dexp=N", "maxC=N", "bounds", "codewords"} <= set(entries(expected)[1])
    for argv in (["--help", "maxC,dest"], ["-h", "maxC,dest"], ["--help=maxC,dest"], ["help=maxC,dest"],
                 ["--help", "maxC", "--morehelp", "dest"]):
        assert dist_m4ri.main(argv) == 0
        assert capsys.readouterr().out == expected, argv


def test_cli_no_args_error(capsys):
    ret = dist_m4ri.main([])
    assert ret == 255
    captured = capsys.readouterr()
    assert "no input matrix or model specified" in captured.err
    assert "Allowed parameters:" in captured.err


def test_cli_version(capsys):
    ret = dist_m4ri.main(["--version"])
    assert ret == 0
    captured = capsys.readouterr()
    assert "0.15.0" in captured.out
    # also the version of the binary
    assert re.search(r"^dist_m4ri version \S+ \(binary '.*dist_m4ri'\)$", captured.out, re.M)


def test_cli_binary_compatibility_silent(capsys):
    ret = dist_m4ri.main(["--help"])
    assert ret == 0
    captured = capsys.readouterr()
    # When binary is found and up to date, stderr should be silent (no warnings)
    assert "Warning:" not in captured.err
    assert "0.15.0" in captured.out


def test_binary_compatibility_warning(tmp_path):
    # Current binary should be compatible and return None
    assert dist_m4ri.check_binary_compatibility() is None

    # Non-existent binary should return a warning
    missing_warn = dist_m4ri.check_binary_compatibility(str(tmp_path / "nonexistent"))
    assert missing_warn is not None
    assert "not found" in missing_warn

    # Older binary script (version 0.5.0) should return a version mismatch warning
    fake_bin = tmp_path / "fake_dist_m4ri"
    fake_bin.write_text("#!/bin/sh\necho \"dist_m4ri version 0.5.0\"\n")
    fake_bin.chmod(0o755)
    older_warn = dist_m4ri.check_binary_compatibility(str(fake_bin))
    assert older_warn is not None
    assert "version 0.5.0" in older_warn
    assert "expected >= 0.15.0" in older_warn


def test_cache_versioning(tmp_path):
    import json
    cache_file = tmp_path / "test_cache.json"

    dist_m4ri.clear_distance_cache()
    dist_m4ri._distance_cache["code_test"] = {"dist": 3, "dmin": 3, "dmax": 3}
    dist_m4ri.save_distance_cache(str(cache_file))

    # Verify version info is silently added to JSON file
    with open(cache_file) as f:
        data = json.load(f)
    assert data.get("__version__") == dist_m4ri.__version__

    # Verify loading does not pollute in-memory code keys
    dist_m4ri.clear_distance_cache()
    dist_m4ri.load_distance_cache(str(cache_file))
    assert "__version__" not in dist_m4ri._distance_cache
    assert "code_test" in dist_m4ri._distance_cache


def test_rw_ksub_kwin_min_hits(capsys):
    dist_m4ri.clear_distance_cache()
    dist_m4ri.disable_distance_cache()
    try:
        dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
        capsys.readouterr()
        d, d_info = dist_m4ri.compute_dem_distance(
            dem=dem_file, method=1, num_steps=2000, ksub=32, kwin=48,
            win_mode=0, min_hits=3, cov_cws=2, refresh=50, threads=4, seed=42
        )
        assert d == 3
        assert d_info[1] == 3
        assert 0 < d_info[2] < 2000
        # the experimental ksub: one warning per call (also when the binary falls back to full-matrix RW)
        assert capsys.readouterr().err.count("# WARNING: ksub=32 is experimental and should not be used") == 1

        args = dist_m4ri.parse_cli_args([
            "method=1", "ksub=64", "win=128", "win_mode=1",
            "min_hits=5", "cov_cws=3", "refresh=100"
        ])
        assert args["ksub"] == 64
        assert args["kwin"] == 128
        assert args["win_mode"] == 1
        assert args["min_hits"] == 5
        assert args["cov_cws"] == 3
        assert args["refresh"] == 100

        default_args = dist_m4ri.parse_cli_args(["method=1"])
        assert default_args["cov_cws"] == 100

        # Test method=3 with min_hits (stops only RW workers while CC proves exact distance)
        d3, d3_info = dist_m4ri.compute_dem_distance(
            dem=dem_file, method=3, num_steps=50000, min_hits=2, cov_cws=2, threads=4, seed=42
        )
        assert d3 == 3
        assert d3_info[0] == 3 and d3_info[1] == 3
        assert "is experimental" not in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()


def test_add_noise_and_noiseless_circuit():
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    # Create a noiseless surface code circuit
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=2,
        distance=3,
        after_clifford_depolarization=0.0
    )
    assert not dist_m4ri.has_noise(circuit)
    noisy = dist_m4ri.add_noise(circuit, p=0.001)
    assert dist_m4ri.has_noise(noisy)

    # compute_dem_distance should automatically add noise if given a noiseless circuit
    dist_m4ri.clear_distance_cache()
    dist, d_info = dist_m4ri.compute_dem_distance(circuit=circuit, threads=4)
    assert dist == 3
    assert d_info[0] == 3 and d_info[1] == 3


def _make_xzzx_circuit(distance: int = 3, rounds: int = 3, memory_basis: str = "Z"):
    """Constructs a rotated-surface-code XZZX circuit by conjugating checkerboard data qubits."""
    import stim
    task = (
        "surface_code:rotated_memory_z"
        if memory_basis.upper() == "Z"
        else "surface_code:rotated_memory_x"
    )
    base_circ = stim.Circuit.generated(
        task,
        rounds=rounds,
        distance=distance,
        after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001,
        after_reset_flip_probability=0.001,
    )
    t_res = dist_m4ri.classify_qubits_thorough(base_circ, basis_arg=memory_basis)
    data_set = set(t_res["data_qubits"])

    coords = {}
    for inst in base_circ.flattened():
        if inst.name == "QUBIT_COORDS":
            args = inst.gate_args_copy()
            for t in inst.targets_copy():
                if t.is_qubit_target:
                    coords[t.value] = args

    rot_data = set()
    for q in data_set:
        x, y = coords[q][0], coords[q][1]
        if int(round((x + y) / 2)) % 2 == 1:
            rot_data.add(q)

    def transform_block(blk):
        out = stim.Circuit()
        for inst in blk:
            if isinstance(inst, stim.CircuitRepeatBlock):
                out.append(
                    stim.CircuitRepeatBlock(
                        inst.repeat_count, transform_block(inst.body_copy())
                    )
                )
            elif inst.name in ["R", "RX"]:
                out.append(inst)
                rot_t = [
                    t.value for t in inst.targets_copy()
                    if t.is_qubit_target and t.value in rot_data
                ]
                if rot_t:
                    out.append("H", rot_t)
            elif inst.name in ["M", "MX"]:
                rot_t = [
                    t.value for t in inst.targets_copy()
                    if t.is_qubit_target and t.value in rot_data
                ]
                if rot_t:
                    out.append("H", rot_t)
                out.append(inst)
            elif inst.name == "CX":
                t_vals = [t.value for t in inst.targets_copy()]
                for i in range(0, len(t_vals), 2):
                    c, t = t_vals[i], t_vals[i + 1]
                    if t in rot_data:
                        out.append("CZ", [c, t])
                    elif c in rot_data:
                        out.append("H", [c])
                        out.append("CX", [c, t])
                        out.append("H", [c])
                    else:
                        out.append("CX", [c, t])
            else:
                out.append(inst)
        return out

    return transform_block(base_circ), rot_data


def test_stim_simple_vs_full_and_basis_tracking(tmp_path, capsys):
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim

    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=3,
        distance=3,
        after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001,
        after_reset_flip_probability=0.001,
    )
    stim_file = tmp_path / "surf_d3_r3_Z.stim"
    circuit.to_file(str(stim_file))

    t_res = dist_m4ri.classify_qubits_thorough(circuit)
    assert t_res["is_css"] is True
    assert t_res["is_rotated_css"] is False
    assert t_res["basis"] == "Z"
    assert len(t_res["data_qubits"]) == 9
    assert len(t_res["x_ancillas"]) == 4
    assert len(t_res["z_ancillas"]) == 4

    simp_circ, stripped, kept = dist_m4ri.strip_minority_detectors(
        circuit, "Z", thorough_res=t_res
    )
    assert stripped == 8
    assert kept == 16
    assert simp_circ.num_detectors == 16

    dist_m4ri.clear_distance_cache()
    dist_m4ri.disable_distance_cache()
    try:
        # The DEM has all 24 detectors; the binary drops the 8 minority-basis detectors (drop=auto, the default)
        ret = dist_m4ri.main([str(stim_file), "-v", "--no-cache", "threads=2"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "Circuit detectors: 24\n" in captured.out  # (the minority-basis detectors only with drop=css)
        assert "dropped 8 of 24 detectors (drop=auto)" in captured.out  # (the binary's summary, printed with -v)
        assert "3 3 0 (exact)" in captured.out
        # drop=css: the minority-basis detectors of the classification, passed to the binary in a file
        ret = dist_m4ri.main([str(stim_file), "drop=css", "-v", "--no-cache", "threads=2"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "Circuit detectors: 24, minority-basis=8" in captured.out
        assert re.search(r"dropped 8 of 24 detectors \(drop=\S+_drop\.txt\)", captured.out)
        assert "3 3 0 (exact)" in captured.out

        # --simple and --full are disabled (since 0.13.0): ignored with a warning, the same DEM
        for opt in ("--simple", "--full", "full=1", "simple=0"):
            ret = dist_m4ri.main([str(stim_file), opt, "-v", "--no-cache", "threads=2"])
            assert ret == 0
            captured = capsys.readouterr()
            name = "--" + opt.split("=")[0].lstrip("-")
            assert f"# Warning: {name}: disabled since version 0.13.0, ignored" in captured.err
            assert "dropped 8 of 24 detectors (drop=auto)" in captured.out
            assert "3 3 0 (exact)" in captured.out
    finally:
        dist_m4ri.enable_distance_cache()


def test_stim_xzzx_rotated_css(tmp_path, capsys):
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")

    for mem_basis in ("Z", "X"):
        xzzx_circ, rot_data = _make_xzzx_circuit(
            distance=3, rounds=3, memory_basis=mem_basis
        )
        t_res = dist_m4ri.classify_qubits_thorough(
            xzzx_circ, basis_arg=mem_basis
        )
        assert t_res["is_css"] is True
        assert t_res["is_rotated_css"] is True
        assert t_res["basis"] == mem_basis
        assert len(t_res["data_qubits"]) == 9
        assert len(t_res["x_ancillas"]) == 4
        assert len(t_res["z_ancillas"]) == 4
        assert set(t_res["rotated_data_qubits"]) == rot_data

        simp_circ, stripped, kept = dist_m4ri.strip_minority_detectors(
            xzzx_circ, mem_basis, thorough_res=t_res
        )
        assert stripped == 8
        assert kept == 16
        assert simp_circ.num_detectors == 16

        stim_path = tmp_path / f"xzzx_d3_r3_{mem_basis}.stim"
        xzzx_circ.to_file(str(stim_path))

        dist_m4ri.clear_distance_cache()
        d_full, info_full = dist_m4ri.compute_dem_distance(dem=str(stim_path), threads=2, debug=1)
        assert d_full == 3
        assert info_full[0] == 3 and info_full[1] == 3
        captured = capsys.readouterr()
        assert "dropped 8 of 24 detectors (drop=auto)" in captured.out  # (the minority-basis detectors)

        # simple= and full= are disabled (since 0.13.0): ignored with a warning to stderr
        for kw in ({"simple": True}, {"full": True}):
            dist_m4ri.clear_distance_cache()
            d, info = dist_m4ri.compute_dem_distance(dem=str(stim_path), threads=2, **kw)
            assert (d, info[0], info[1]) == (3, 3, 3)
            assert f"# Warning: {next(iter(kw))}=: disabled since version 0.13.0" in capsys.readouterr().err


def test_stim_rounds_option(tmp_path, capsys):
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim

    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=5,
        distance=3,
        after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001,
        after_reset_flip_probability=0.001,
    )
    stim_with_repeat = tmp_path / "surf_r5_Z.stim"
    circuit.to_file(str(stim_with_repeat))

    # CLI argument parsing checks for --rounds, --simple, --full (the last two are disabled, parsed for the warning)
    args1 = dist_m4ri.parse_cli_args(["--rounds", str(stim_with_repeat)])
    assert args1["rounds"] == 2
    assert args1["fdem"] == str(stim_with_repeat)

    args2 = dist_m4ri.parse_cli_args(["--rounds", "4", "--simple", str(stim_with_repeat)])
    assert args2["rounds"] == 4
    assert args2["simple"] is True
    assert args2["full"] is None
    assert args2["fdem"] == str(stim_with_repeat)

    args3 = dist_m4ri.parse_cli_args(["rounds=3", "--full", "basis=X", str(stim_with_repeat)])
    assert args3["rounds"] == 3
    assert args3["full"] is True
    assert args3["simple"] is None
    assert args3["basis"] == "X"

    dist_m4ri.clear_distance_cache()
    dist_m4ri.disable_distance_cache()
    try:
        # Bare --rounds sets rounds=2 (REPEAT 2 block -> 24 detectors, of which 8 X-detectors are dropped)
        ret = dist_m4ri.main(
            [str(stim_with_repeat), "--rounds", "-v", "--no-cache", "threads=2"]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "rounds=2" in out
        assert "Circuit detectors: 24\n" in out
        assert "dropped 8 of 24 detectors (drop=auto)" in out

        # Flattened circuit has no REPEAT block -> issues a warning and continues
        flat_stim = tmp_path / "surf_flat_Z.stim"
        circuit.flattened().to_file(str(flat_stim))
        d_flat, info_flat = dist_m4ri.compute_dem_distance(
            dem=str(flat_stim), rounds=2, threads=2
        )
        assert d_flat == 3
        assert info_flat[0] == 3 and info_flat[1] == 3
        err_api = capsys.readouterr().err
        assert "Warning:" in err_api and "REPEAT block" in err_api

        ret_flat = dist_m4ri.main(
            [str(flat_stim), "--rounds", "--no-cache", "threads=2"]
        )
        assert ret_flat == 0
        captured_flat = capsys.readouterr()
        assert "Warning:" in captured_flat.err and "REPEAT block" in captured_flat.err

        # Non-stim input with --rounds must signal an error
        dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
        ret_dem_err = dist_m4ri.main([dem_file, "--rounds", "--no-cache"])
        assert ret_dem_err != 0
        assert "--rounds" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()


def test_stim_out_dem_out_stim_out_dir(tmp_path):
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim

    in_dir = tmp_path / "inputs"
    in_dir.mkdir(parents=True, exist_ok=True)
    stim_file = in_dir / "surf_d3_Z.stim"

    # Create a noiseless circuit to also verify that --out-stim saves the noisy/processed circuit
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=4,
        distance=3,
    )
    circuit.to_file(str(stim_file))

    # 1. Test CLI parsing of bare --out-dem and --out-stim vs explicit filenames
    p1 = dist_m4ri.parse_cli_args(["--out-dem", "--out-stim", str(stim_file)])
    assert p1["out_dem"] is True
    assert p1["out_stim"] is True
    assert p1["fdem"] == str(stim_file)

    p2 = dist_m4ri.parse_cli_args(
        ["--out-dir", "/tmp/out", "--out-dem", "a.dem", "--out-stim", "b.stim", str(stim_file)]
    )
    assert p2["out_dir"] == "/tmp/out"
    assert p2["out_dem"] == "a.dem"
    assert p2["out_stim"] == "b.stim"
    assert p2["fdem"] == str(stim_file)

    dist_m4ri.clear_distance_cache()
    dist_m4ri.disable_distance_cache()
    try:
        # 2. Bare --out-dem and --out-stim (all detectors, default out-dir = input file dir)
        ret = dist_m4ri.main(
            [str(stim_file), "--rounds", "2", "--out-dem", "--out-stim", "--no-cache", "threads=2"]
        )
        assert ret == 0
        exp_dem = in_dir / "surf_d3_Z_full.dem"
        exp_stim = in_dir / "surf_d3_Z_full.stim"
        assert exp_dem.exists()
        assert exp_stim.exists()

        saved_circ = stim.Circuit.from_file(str(exp_stim))
        assert dist_m4ri.has_noise(saved_circ)
        assert saved_circ.num_detectors == 24
        saved_dem = stim.DetectorErrorModel.from_file(str(exp_dem))
        assert saved_dem.num_detectors == 24

        # 3. Bare --out-dem and --out-stim with explicit --out-dir
        custom_dir = tmp_path / "custom_out"
        ret_full = dist_m4ri.main(
            [
                str(stim_file),
                "--rounds",
                "2",
                "--out-dir",
                str(custom_dir),
                "--out-dem",
                "--out-stim",
                "--no-cache",
                "threads=2",
            ]
        )
        assert ret_full == 0
        exp_full_dem = custom_dir / "surf_d3_Z_full.dem"
        exp_full_stim = custom_dir / "surf_d3_Z_full.stim"
        assert exp_full_dem.exists()
        assert exp_full_stim.exists()
        saved_full_dem = stim.DetectorErrorModel.from_file(str(exp_full_dem))
        assert saved_full_dem.num_detectors == 24

        # 4. Explicit bare filenames without --out-dir default to input file's directory
        ret_named_in_dir = dist_m4ri.main(
            [
                str(stim_file),
                "--rounds",
                "2",
                "--out-dem",
                "named_in_dir.dem",
                "--out-stim",
                "named_in_dir.stim",
                "--no-cache",
                "threads=2",
            ]
        )
        assert ret_named_in_dir == 0
        assert (in_dir / "named_in_dir.dem").exists()
        assert (in_dir / "named_in_dir.stim").exists()

        # 5. Explicit bare filenames with --out-dir save into --out-dir
        ret_named = dist_m4ri.main(
            [
                str(stim_file),
                "--rounds",
                "2",
                "--out-dir",
                str(custom_dir),
                "--out-dem",
                "explicit.dem",
                "--out-stim",
                "explicit.stim",
                "--no-cache",
                "threads=2",
            ]
        )
        assert ret_named == 0
        assert (custom_dir / "explicit.dem").exists()
        assert (custom_dir / "explicit.stim").exists()
    finally:
        dist_m4ri.enable_distance_cache()


S5_H = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
S5_L = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")


def test_normalize_start_list():
    norm = dist_m4ri._normalize_start_list
    assert norm(None) is None
    assert norm(5) == [5]
    assert norm(-1) is None
    assert norm("-1") is None
    assert norm("") is None
    assert norm("48,0,48") == [0, 48]
    assert norm([7, 3, 3]) == [3, 7]
    assert norm((2,)) == [2]
    assert norm(np.array([4, 1])) == [1, 4]
    for bad in ("1,,2", "a", "1.5", "0,-1", [1, -2], 1.5, True, [True]):
        with pytest.raises(ValueError):
            norm(bad)
    assert dist_m4ri._start_key_suffix([0, 48]) == ":start=0,48"
    assert dist_m4ri._start_key_suffix(None) == ""


def test_run_dist_m4ri_start_list_and_disabled_options(capsys):
    # start list is passed sorted and deduplicated; noscan/cbeg/cend are ignored with a warning
    res = dist_m4ri.run_dist_m4ri(
        method=2, finH=S5_H, finL=S5_L, wmax=7, start=[41, 0, 41], noscan=1, cbeg=3, cend=7,
        threads=2, verbose=True
    )
    assert res == (5, 5, 0)
    captured = capsys.readouterr()
    run_line = [ln for ln in captured.out.splitlines() if "[dist_m4ri] Running:" in ln][0]
    assert "start=0,41" in run_line
    assert "noscan=" not in run_line and "cbeg=" not in run_line and "cend=" not in run_line
    for opt in ("noscan=1", "cbeg=3", "cend=7"):
        assert f"{opt} is disabled in the Python interface" in captured.err
    assert "WARNING: start=0,41 (expert option)" in captured.err


def test_sort_rows_option(tmp_path, capsys):
    # sort_rows=None leaves the binary's default (rows of H sorted by weight for matrices, not for a DEM); True / False
    # pass sort_rows=1 / sort_rows=0; the distance and the codewords enumerated (outC) do not depend on the row order
    cws = {}
    for val, opt in ((None, None), (True, "sort_rows=1"), (False, "sort_rows=0")):
        out = str(tmp_path / f"cws_{opt}.nz")
        res = dist_m4ri.run_dist_m4ri(method=2, finH=S5_H, finL=S5_L, wmax=5, threads=2, outC=out, sort_rows=val,
                                      verbose=True)
        assert tuple(res) == (5, 5, 0)
        captured = capsys.readouterr()
        run_line = [ln for ln in captured.out.splitlines() if "[dist_m4ri] Running:" in ln][0]
        if opt is None:
            assert "sort_rows=" not in run_line
        else:
            assert opt in run_line
        sorted_note = "rows of H sorted by increasing weight for CC" in captured.out + captured.err
        assert sorted_note == (val is not False)
        cws[opt] = sorted(tuple(cw) for cw in dist_m4ri.read_sparse_vectors(out))
    assert len(cws[None]) > 1 and cws[None] == cws["sort_rows=1"] == cws["sort_rows=0"]
    # the compute_*_distance() functions and the CLI pass it on (the cache key does not depend on it)
    dist_m4ri.clear_distance_cache()
    try:
        for val in (True, False):
            assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, method=2, wmax=5, threads=2, sort_rows=val,
                                                      debug=dist_m4ri.PY_DBG_COMMANDS, cache_file=None) == 5
            assert f"sort_rows={int(val)}" in capsys.readouterr().out
            dist_m4ri.clear_distance_cache()
        assert dist_m4ri.parse_cli_args(["sort_rows=0"])["sort_rows"] is False
        assert dist_m4ri.parse_cli_args(["sort_rows=1"])["sort_rows"] is True
        assert dist_m4ri.parse_cli_args([])["sort_rows"] is None
        assert dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "method=2", "wmax=5", "threads=2", "sort_rows=0",
                               "--no-cache", f"debug={dist_m4ri.PY_DBG_COMMANDS}"]) == 0
        assert "sort_rows=0" in capsys.readouterr().out
        assert dist_m4ri.main([f"finH={S5_H}", "sort_rows=2", "--no-cache"]) == 1
        assert "invalid 'sort_rows=2': sort_rows should be 0 or 1" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()
        dist_m4ri.clear_distance_cache()


def test_cc_prune_option(tmp_path, capsys):
    # cc_prune=None leaves the binary's default (1, the branch and bound); 0 / 1 / 2 are passed as cc_prune=N; the
    # distance and the codewords enumerated (outC) do not depend on it, and with cc_prune=1 or 2 each codeword is found
    # exactly once (surf_d5 as matrices: with cc_prune=0, codewords with equal columns are found several times)
    cws = {}
    for val in (None, 0, 1, 2):
        out = str(tmp_path / f"cws_{val}.nz")
        res = dist_m4ri.run_dist_m4ri(method=2, finH=S5_H, finL=S5_L, wmax=5, threads=2, outC=out, cc_prune=val,
                                      verbose=True)
        assert tuple(res) == (5, 5, 0)
        captured = capsys.readouterr()
        run_line = [ln for ln in captured.out.splitlines() if "[dist_m4ri] Running:" in ln][0]
        if val is None:
            assert "cc_prune=" not in run_line
        else:
            assert f"cc_prune={val}" in run_line
        hits = re.search(r"# codewords accumulated: total=(\d+), min_w=5: cws=\d+, total_hits=(\d+),",
                         captured.out + captured.err)
        assert hits, captured.err
        assert (hits.group(1) == hits.group(2)) == (val != 0)
        cws[val] = sorted(tuple(cw) for cw in dist_m4ri.read_sparse_vectors(out))
    assert len(cws[None]) > 1 and cws[None] == cws[0] == cws[1] == cws[2]
    with pytest.raises(ValueError, match="cc_prune should be 0, 1, or 2"):
        dist_m4ri.run_dist_m4ri(method=2, finH=S5_H, finL=S5_L, wmax=5, cc_prune=3)
    # the compute_*_distance() functions and the CLI pass it on (the cache key does not depend on it)
    dist_m4ri.clear_distance_cache()
    try:
        for val in (0, 2):
            assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, method=2, wmax=5, threads=2, cc_prune=val,
                                                      debug=dist_m4ri.PY_DBG_COMMANDS, cache_file=None) == 5
            assert f"cc_prune={val}" in capsys.readouterr().out
            dist_m4ri.clear_distance_cache()
        assert dist_m4ri.parse_cli_args(["cc_prune=0"])["cc_prune"] == 0
        assert dist_m4ri.parse_cli_args(["--cc-prune", "2"])["cc_prune"] == 2
        assert dist_m4ri.parse_cli_args([])["cc_prune"] is None
        assert dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "method=2", "wmax=5", "threads=2", "--cc-prune", "0",
                               "--no-cache", f"debug={dist_m4ri.PY_DBG_COMMANDS}"]) == 0
        assert "cc_prune=0" in capsys.readouterr().out
        assert dist_m4ri.main([f"finH={S5_H}", "cc_prune=3", "--no-cache"]) == 1
        assert "invalid 'cc_prune=3': cc_prune should be 0, 1, or 2" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()
        dist_m4ri.clear_distance_cache()


def test_cc_round_stats(capsys):
    # with debug=8 (bit 8 of the binary), _last_run_stats["cc_rounds"] lists the completed CC rounds (weight, CC thread
    # time, predicted work, wall time, CC threads, and nodes: the calls of the CC recursion), and "cc_nodes" is the
    # total of the nodes; the lines of binaries before 0.15.0 have no nodes (no "cc_nodes").  The number of nodes of a
    # completed round does not depend on the threads (surf_d5: the rounds w=1..4 are completed, the round w=5 stops at
    # the first codeword; method=2 or method=3 without RW steps)
    stats = dist_m4ri._parse_stderr_stats(
        "# CC round w=3: started 2 more threads, 4 in total (work not known)\n"
        "# CC round w=3 work: 0.00125 thread-s measured, unknown predicted (0.001s on 4 threads), 1956 nodes\n"
        "# CC round w=4 started: 8 CC threads, 0 RW threads (bounds [4, 0], rem_rw=0, rem_time=0.00s, est 0.002 "
        "thread-s)\n"
        "# CC round w=4 work: 1.5e-05 thread-s measured, 0.002 thread-s predicted (0.004s, 8 CC threads at the end), "
        "24644 nodes\n")
    assert stats["cc_rounds"] == [
        {"w": 3, "work": 0.00125, "predicted": None, "time": 0.001, "threads": 4, "nodes": 1956},
        {"w": 4, "work": 1.5e-05, "predicted": 0.002, "time": 0.004, "threads": 8, "nodes": 24644}]
    assert stats["cc_nodes"] == 1956 + 24644
    old = dist_m4ri._parse_stderr_stats(
        "# CC round w=3 work: 0.001 thread-s measured, unknown predicted (0.001s on 2 threads)\n")
    assert old["cc_rounds"][0]["nodes"] is None and "cc_nodes" not in old
    assert "cc_rounds" not in dist_m4ri._parse_stderr_stats("# codewords accumulated: total=0\n")
    nodes = {}
    for method, threads in ((2, 1), (2, 2), (3, 2)):
        res = dist_m4ri.run_dist_m4ri(method=method, finH=S5_H, finL=S5_L, wmax=5, threads=threads, debug=8,
                                      steps=0 if method == 3 else None)
        assert tuple(res) == (5, 5, 0)
        capsys.readouterr()
        stats = dist_m4ri._last_run_stats
        rounds = stats["cc_rounds"]
        assert [r["w"] for r in rounds] == [1, 2, 3, 4]
        assert all(1 <= r["threads"] <= threads and r["time"] >= 0.0 and r["work"] >= 0.0 for r in rounds)
        assert rounds[0]["nodes"] == 0 and all(r["nodes"] > 0 for r in rounds[1:])
        assert stats["cc_nodes"] == sum(r["nodes"] for r in rounds)
        nodes[(method, threads)] = [r["nodes"] for r in rounds]
    assert nodes[(2, 1)] == nodes[(2, 2)] == nodes[(3, 2)]
    dist_m4ri.run_dist_m4ri(method=2, finH=S5_H, finL=S5_L, wmax=5, threads=2, debug=1)
    capsys.readouterr()
    assert "cc_rounds" not in dist_m4ri._last_run_stats and "cc_nodes" not in dist_m4ri._last_run_stats


def test_win_hyb_option(capsys):
    # win_hyb (window size of the hybrid RW steps, n >= 500; 0, the default: 256 columns for a DEM, 512 for matrices)
    # is passed to the binary if nonzero by run_dist_m4ri(), the compute_*_distance() functions, and the CLI; the
    # input summary (debug&1) reports the windows, and the codewords found by the uniform and by the window steps are
    # parsed (_last_run_stats["rw_step_types"]), as the rounds of a DEM (_last_run_stats["rounds"])
    for val, note in ((0, "512 columns on every other step, uniform permutations on the others (win_hyb=0: 256 for a "
                          "DEM, 512 for matrices)"),
                      (100, "100 columns on every other step, uniform permutations on the others (win_hyb=100)")):
        res = dist_m4ri.run_dist_m4ri(method=1, finH=S5_H, finL=S5_L, steps=200, threads=2, win_hyb=val, verbose=True)
        assert res[2] > 1
        captured = capsys.readouterr()
        run_line = [ln for ln in captured.out.splitlines() if "[dist_m4ri] Running:" in ln][0]
        assert ("win_hyb=" in run_line) == (val != 0)
        windows = [ln for ln in (captured.out + captured.err).splitlines() if "RW windows:" in ln]
        assert len(windows) == 1 and windows[0].endswith(note)
        types = dist_m4ri._last_run_stats["rw_step_types"]
        assert types and all(t["uniform_cws"] + t["window_cws"] - t["both_cws"] == t["cws"] for t in types.values())
        assert "rounds" not in dist_m4ri._last_run_stats  # (a matrix)
        assert "by RW step type" in dist_m4ri.explain_bounds(list(res), method=1, stats=dist_m4ri._last_run_stats)
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    dist_m4ri.run_dist_m4ri(method=2, fdem=dem_file, wmax=3, threads=2, debug=1)
    capsys.readouterr()
    rounds = dist_m4ri._last_run_stats["rounds"]
    assert set(rounds) == {"D", "L", "score"} and 0.0 <= rounds["score"] <= 1.0
    dist_m4ri.clear_distance_cache()
    try:
        assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, method=1, num_steps=20, min_hits=0, threads=2,
                                                  win_hyb=300, debug=dist_m4ri.PY_DBG_COMMANDS, cache_file=None) > 0
        assert "win_hyb=300" in capsys.readouterr().out
        dist_m4ri.clear_distance_cache()
        assert dist_m4ri.parse_cli_args(["win_hyb=200"])["win_hyb"] == 200
        assert dist_m4ri.parse_cli_args([])["win_hyb"] == 0
        assert dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "method=1", "steps=20", "threads=2", "win_hyb=100",
                               "--no-cache", f"debug={dist_m4ri.PY_DBG_COMMANDS}"]) == 0
        assert "win_hyb=100" in capsys.readouterr().out
        assert dist_m4ri.main([f"finH={S5_H}", "method=1", "win_hyb=-1", "--no-cache"]) == 1
        assert "invalid 'win_hyb=-1'" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()
        dist_m4ri.clear_distance_cache()


def test_pairs_option(capsys):
    # pairs (the RW pairs of non-pivot columns, p=2: the pairs of columns which agree on N random pivot rows; -1, the
    # default: auto; 0: off) is passed to the binary if not -1 by run_dist_m4ri(), the compute_*_distance() functions,
    # and the CLI; the summary (debug&1) gives the number l of rows (auto: 8 for surf_d5 with n=1958, rank 120, and
    # d=5), the pairs checked per step, and the codewords found by pairs and by single columns (_last_run_stats)
    for val, l_exp in ((-1, 8), (0, 0), (12, 12)):
        res = dist_m4ri.run_dist_m4ri(method=1, finH=S5_H, finL=S5_L, steps=200, threads=1, seed=1, pairs=val,
                                      verbose=True)
        assert res[1] == 5
        captured = capsys.readouterr()
        run_line = [ln for ln in captured.out.splitlines() if "[dist_m4ri] Running:" in ln][0]
        assert ("pairs=" in run_line) == (val != -1)
        stats = dist_m4ri._last_run_stats
        assert stats["rw_pairs_l"] == l_exp
        assert ("rw_pairs_per_step" in stats) == (l_exp > 0)
        exp = dist_m4ri.explain_bounds(list(res), method=1, stats=stats)
        assert (f"RW pairs (p=2, buckets of l = {l_exp} pivot rows): " in exp) == (l_exp > 0)
    dist_m4ri.clear_distance_cache()
    try:
        assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, method=1, num_steps=20, min_hits=0, threads=2,
                                                  pairs=0, debug=dist_m4ri.PY_DBG_COMMANDS, cache_file=None) > 0
        assert "pairs=0" in capsys.readouterr().out
        dist_m4ri.clear_distance_cache()
        assert dist_m4ri.parse_cli_args(["pairs=4"])["pairs"] == 4
        assert dist_m4ri.parse_cli_args([])["pairs"] == -1
        assert dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "method=1", "steps=20", "threads=2", "pairs=6",
                               "--no-cache", f"debug={dist_m4ri.PY_DBG_COMMANDS}"]) == 0
        assert "pairs=6" in capsys.readouterr().out
        for bad in (-2, 65):
            assert dist_m4ri.main([f"finH={S5_H}", "method=1", f"pairs={bad}", "--no-cache"]) == 1
            assert f"invalid 'pairs={bad}'" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()
        dist_m4ri.clear_distance_cache()


def test_model_lower_bound(tmp_path, capsys):
    # the code-sized model of a Stim memory circuit (compute_model_lower_bound(), CLI --model-lb): one row per chain of
    # detectors linked by their measurements (the checks over the rounds; by default only the chains with final data
    # measurements, all of them with minority=True, CLI --model-minority), the sum of their rows, so that every DEM
    # codeword gives a model codeword of no larger weight and the model distance is a certified lower bound: 3 (tight)
    # for the surface code with d=3 (8 checks, 4 of the measured basis), 5 for the bb_72 circuit of distance 6 (also
    # with the minority-basis rows); the bound is stored as dmin in the cache record
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circ = stim.Circuit.generated("surface_code:rotated_memory_x", distance=3, rounds=3,
                                  after_clifford_depolarization=0.001, before_measure_flip_probability=0.001,
                                  after_reset_flip_probability=0.001, before_round_data_depolarization=0.001)
    path = str(tmp_path / "surf_x_d3.stim")
    circ.to_file(path)
    cache = str(tmp_path / "cache.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()  # (a CLI test with --no-cache may have disabled it)
    try:
        b, info = dist_m4ri.compute_model_lower_bound(path, threads=2, cache_file=cache, out_model=str(tmp_path / "m"),
                                                      return_info=True)
        assert b == 3 and info["model_d_info"][:2] == [3, 3] and info["d_info"] == [3, 0, 0]
        assert (info["groups"], info["groups_final"], info["rows"], info["rows_measured"]) == (8, 4, 4, 4)
        assert info["trivial"] == 0 and len(info["col_map"]) == info["dem_cols"]
        assert info["files"] == [str(tmp_path / "m_HR.mtx"), str(tmp_path / "m_LR.mtx")]
        assert all(os.path.isfile(f) for f in info["files"])
        entry = dist_m4ri.get_cached_distance(circuit=path, cache_file=cache)
        assert (entry["dmin"], entry["dmax"], entry["dmin_model"]) == (3, 0, 3)
        b, info = dist_m4ri.compute_model_lower_bound(path, threads=2, cache_file=cache, minority=True,
                                                      return_info=True)
        assert b == 3 and (info["rows"], info["rows_measured"]) == (8, 4)
        # the DEM run starts CC from the cached bound; its codewords map to model codewords of no larger weight
        dist, d_info, cws = dist_m4ri.compute_dem_distance(circuit=path, method=2, wmax=4, dW=1, drop=0, do_cws=True,
                                                           threads=2, cache_file=cache)
        assert dist == 3 and d_info[:2] == [3, 3] and len(cws) > 0
        entry = dist_m4ri.get_cached_distance(circuit=path, cache_file=cache)
        assert (entry["dmin"], entry["dmax"], entry["dmin_model"]) == (3, 3, 3)
        H, L = info["H"].toarray(), info["L"].toarray()
        for cw in cws:
            y = np.zeros(H.shape[1], dtype=int)
            for j in cw:
                if info["col_map"][j] >= 0:
                    y[info["col_map"][j]] ^= 1
            assert not (H @ y % 2).any() and (L @ y % 2).any() and y.sum() <= len(cw)
        # a circuit in the repository: the model is 1 below the distance 6 (also with the measured-basis rows only)
        bb72 = os.path.join(EXAMPLES_DIR, "..", "benchmark", "stim_bb", "bb_72_12_6_si1000_r6_X.stim")
        for minority in (True, False):
            b, info = dist_m4ri.compute_model_lower_bound(bb72, threads=2, minority=minority, cache_file=None,
                                                          return_info=True)
            assert b == 5 and info["rows"] == (72 if minority else 36) and info["rows_measured"] == 36
        # an error which flips the observable but no detector: the bound is 1, with a note
        triv = str(tmp_path / "triv.stim")
        with open(triv, "w") as f:
            f.write("R 0 1 2 3\nX_ERROR(0.1) 0 1 2\nCX 0 3 1 3\nMR 3\nDETECTOR rec[-1]\nM 0 1 2\n"
                    "DETECTOR rec[-2] rec[-3] rec[-4]\nOBSERVABLE_INCLUDE(0) rec[-1]\n")
        b, info = dist_m4ri.compute_model_lower_bound(triv, threads=2, cache_file=None, return_info=True)
        assert b == 1 and info["trivial"] == 1 and (info["groups"], info["groups_final"]) == (1, 1)
        assert "flip an observable but no group of detectors" in capsys.readouterr().err
        # sort_rows and cc_prune are passed to each run of the binary (the bound does not depend on them)
        dist_m4ri.clear_distance_cache()  # (the model of the first runs is cached)
        b = dist_m4ri.compute_model_lower_bound(path, threads=2, cache_file=None, minority=True, sort_rows=False,
                                                cc_prune=0, debug=dist_m4ri.PY_DBG_COMMANDS)
        runs = [ln for ln in capsys.readouterr().out.splitlines() if "[dist_m4ri] Running:" in ln]
        assert b == 3 and len(runs) >= 2 and all("sort_rows=0" in ln and "cc_prune=0" in ln for ln in runs)
        # CLI: the bounds of the circuit (without the cache: the bound alone), the model files; errors and warnings
        assert dist_m4ri.main([f"fdem={path}", "--model-lb", "--no-cache", "threads=2", "--out-model",
                               str(tmp_path / "cli")]) == 0
        assert capsys.readouterr().out.strip().splitlines()[-1] == "3 0 0"
        assert os.path.isfile(str(tmp_path / "cli_HR.mtx")) and os.path.isfile(str(tmp_path / "cli_LR.mtx"))
        from scipy.io import mmread
        assert mmread(str(tmp_path / "cli_HR.mtx")).shape[0] == 4  # (the measured-basis rows only, the default)
        assert dist_m4ri.main([f"fdem={path}", "--model-lb", "--model-minority", "--no-cache", "threads=2",
                               "--out-model", str(tmp_path / "clim")]) == 0
        assert capsys.readouterr().out.strip().splitlines()[-1] == "3 0 0"
        assert mmread(str(tmp_path / "clim_HR.mtx")).shape[0] == 8  # (also the minority-basis rows)
        assert dist_m4ri.main([f"fdem={path}", "--model-lb", "--no-cache", "threads=2", "sort_rows=1", "--cc-prune",
                               "2", f"debug={dist_m4ri.PY_DBG_COMMANDS}"]) == 0
        out = capsys.readouterr().out.strip().splitlines()
        runs = [ln for ln in out if "[dist_m4ri] Running:" in ln]
        assert out[-1] == "3 0 0" and runs and all("sort_rows=1" in ln and "cc_prune=2" in ln for ln in runs)
        assert dist_m4ri.main([f"fdem={path}", "model_minority=1", "--no-cache", "threads=2", "method=2"]) == 0
        assert "# Warning: --model-minority: only used with --model-lb; ignored" in capsys.readouterr().err
        dist_m4ri.enable_distance_cache()  # (--no-cache disabled it)
        assert dist_m4ri.main([f"fdem={path}", "--model-lb", f"cache={cache}", "threads=2", "outC=x.nz"]) == 0
        captured = capsys.readouterr()
        assert captured.out.strip().splitlines()[-1] == "3 3 0 (exact)"
        assert "# Warning: outC=: not used with --model-lb; ignored" in captured.err
        dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
        assert dist_m4ri.main([f"fdem={dem_file}", "--model-lb", "--no-cache"]) == 1
        assert "--model-lb needs a Stim circuit" in capsys.readouterr().err
        with pytest.raises(ValueError, match="not a Stim circuit"):
            dist_m4ri.compute_model_lower_bound(dem_file)
        with pytest.raises(TypeError, match="did you mean 'minority'"):
            dist_m4ri.compute_model_lower_bound(path, minorty=False)
    finally:
        dist_m4ri.enable_distance_cache()
        dist_m4ri.clear_distance_cache()


def test_start_list_cache_records(tmp_path, monkeypatch, capsys):
    import json
    cache_file = str(tmp_path / "start_cache.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    try:
        # 1. Start-list run: separate ':start=0,41' record; the main record only gets dmax (and codewords)
        dist, d_info = dist_m4ri.compute_quantum_distance(
            S5_H, L=S5_L, method=2, wmax=7, start="41,0", threads=2, return_info=True, cache_file=cache_file
        )
        assert dist == 5 and d_info == [5, 5, 0]
        err = capsys.readouterr().err
        assert "WARNING: start=0,41 (expert option)" in err
        assert "use trust_start=1 to also accept dmin" in err
        plain = dist_m4ri.get_cached_distance(H=S5_H, L=S5_L, cache_file=cache_file)
        srec = dist_m4ri.get_cached_distance(H=S5_H, L=S5_L, cache_file=cache_file, start=[0, 41])
        assert srec["dmin"] == 5 and srec["dmax"] == 5
        assert plain["dmin"] == 0 and plain["dmax"] == 5
        with open(cache_file) as f:
            data = json.load(f)
        assert any(k.startswith("quantum:") and k.endswith(":start=0,41") for k in data)

        # 2. Repeated start-list call is a cache hit (the binary is not called)
        def _no_run(*args, **kwargs):
            raise AssertionError("dist_m4ri binary should not be called")
        monkeypatch.setattr(dist_m4ri, "run_dist_m4ri", _no_run)
        dist2, d_info2 = dist_m4ri.compute_quantum_distance(
            S5_H, L=S5_L, method=2, wmax=7, start=[0, 41], threads=2, return_info=True, cache_file=cache_file
        )
        assert dist2 == 5 and d_info2 == [5, 5, 0]
        assert dist_m4ri.get_cached_distance(H=S5_H, L=S5_L, cache_file=cache_file)["dmin"] == 0

        # 3. trust_start copies the exact start-list record to the main record without a calculation
        dist3, d_info3 = dist_m4ri.compute_quantum_distance(
            S5_H, L=S5_L, method=2, wmax=7, start=[0, 41], trust_start=True, threads=2, return_info=True,
            cache_file=cache_file
        )
        assert dist3 == 5 and d_info3 == [5, 5, 0]
        assert "trust_start=1: start-list results are accepted" in capsys.readouterr().err
        with open(cache_file) as f:
            data = json.load(f)
        plain_keys = [k for k in data if k.startswith("quantum:") and ":start=" not in k]
        assert len(plain_keys) == 1
        assert data[plain_keys[0]]["dmin"] == 5 and data[plain_keys[0]]["dmax"] == 5

        # 4. The main record is now exact: a call without a start list is a cache hit
        dist4, d_info4 = dist_m4ri.compute_quantum_distance(
            S5_H, L=S5_L, method=2, wmax=7, threads=2, return_info=True, cache_file=cache_file
        )
        assert dist4 == 5 and d_info4 == [5, 5, 0]
    finally:
        dist_m4ri.clear_distance_cache()


def test_start_list_trust_after_calculation_css(tmp_path):
    cache_file = str(tmp_path / "start_css_cache.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    try:
        kw = dict(Hx=S5_H, Hz=S5_H, Lz=S5_L, Lx=S5_L, method=2, wmax=7, threads=2, cache_file=cache_file)
        # Without trust_start: only the upper bounds (and codewords) reach the main record
        dist, dx, dz = dist_m4ri.compute_css_distance(start="0,41", **kw)
        assert dist == 5 and dx == [5, 5, 0] and dz == [5, 5, 0]
        plain = dist_m4ri.get_cached_distance(Hx=S5_H, Hz=S5_H, Lx=S5_L, Lz=S5_L, cache_file=cache_file)
        assert plain["dmax_X"] == 5 and plain["dmax_Z"] == 5
        assert plain["dmin_X"] == 0 and plain["dmin_Z"] == 0
        # With trust_start and another start list: after the calculation, dmin is merged as well
        dist2, dx2, dz2 = dist_m4ri.compute_css_distance(start=[0, 41, 12], trust_start=True, **kw)
        assert dist2 == 5 and dx2 == [5, 5, 0] and dz2 == [5, 5, 0]
        plain2 = dist_m4ri.get_cached_distance(Hx=S5_H, Hz=S5_H, Lx=S5_L, Lz=S5_L, cache_file=cache_file)
        assert plain2["dX"] == [5, 5, 0] and plain2["dZ"] == [5, 5, 0]
        assert plain2["dmin"] == 5 and plain2["dmax"] == 5 and plain2["dist"] == 5
    finally:
        dist_m4ri.clear_distance_cache()


def test_cli_start_trust_start_and_disabled_options(capsys):
    args = dist_m4ri.parse_cli_args(
        ["finH=h.mtx", "start=5,1,1", "--trust-start", "noscan=1", "cbeg=2", "cend=3", "irred=0"]
    )
    assert args["start"] == [1, 5]
    assert args["trust_start"] is True
    assert (args["noscan"], args["cbeg"], args["cend"], args["irred"]) == (1, 2, 3, 0)
    assert dist_m4ri.parse_cli_args([])["irred"] is None
    assert dist_m4ri.parse_cli_args(["trust_start=0"])["trust_start"] is False
    assert dist_m4ri.parse_cli_args(["trust-start=1"])["trust_start"] is True
    assert dist_m4ri.parse_cli_args(["--trust-start=yes"])["trust_start"] is True
    assert dist_m4ri.parse_cli_args(["start=-1"])["start"] is None
    assert dist_m4ri.parse_cli_args([])["trust_start"] is False

    try:
        # Invalid start list is an error
        ret = dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "start=1,,2", "--no-cache"])
        assert ret == 1
        assert "invalid start='1,,2'" in capsys.readouterr().err

        # Disabled options are ignored with a warning; the start list gives the distance with a warning
        ret = dist_m4ri.main([
            f"finH={S5_H}", f"finL={S5_L}", "method=2", "wmax=7", "noscan=1", "cbeg=30",
            "start=0,41", "--no-cache", "threads=2", "irred=0"
        ])
        assert ret == 0
        captured = capsys.readouterr()
        assert "5 5 0 (exact)" in captured.out
        assert "noscan=1 is disabled in the Python interface" in captured.err
        assert "cbeg=30 is disabled in the Python interface" in captured.err
        assert "irred=0 is disabled in the Python interface" in captured.err
        assert "WARNING: start=0,41 (expert option)" in captured.err
        # irred=1 (the binary's default) is accepted silently
        ret = dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "method=2", "wmax=7", "--no-cache", "irred=1"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "5 5 0 (exact)" in captured.out
        assert "irred" not in captured.err
    finally:
        dist_m4ri.enable_distance_cache()


TRY_X = os.path.join(EXAMPLES_DIR, "tryX.mtx")  # small CSS code [[40,10,4]]
TRY_Z = os.path.join(EXAMPLES_DIR, "tryZ.mtx")


def test_css_distance_requires_both_matrices(capsys):
    empty = np.zeros((0, 40), dtype=np.int8)
    for hx, hz in ((TRY_X, None), (None, TRY_Z), (TRY_X, empty), (None, None)):
        with pytest.raises(ValueError, match="requires both Hx and Hz"):
            dist_m4ri.compute_css_distance(Hx=hx, Hz=hz, method=2, wmax=4, threads=2)
    try:
        # CLI: an error message and exit code 1
        ret = dist_m4ri.main([f"Hx={TRY_X}", "method=2", "wmax=4", "--no-cache"])
        assert ret == 1
        assert "requires both Hx and Hz" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()


def test_distance_result_eq():
    r = dist_m4ri.DistanceResult(5, 5, 0)
    assert r == 5 and r == np.int64(5) and r != 4
    assert r == (5, 5, 0) and r == [5, 5, 0]
    assert r == dist_m4ri.DistanceResult(5, 5, 0) and r != dist_m4ri.DistanceResult(5, 5, 10)
    assert r != "5"
    assert len(r) == 3 and r[1] == 5 and tuple(r) == (5, 5, 0)
    rc = dist_m4ri.DistanceResult(4, 4, 0, cws=[[0, 1, 2, 3]])
    assert len(rc) == 4 and rc[3] == [[0, 1, 2, 3]] and rc == (4, 4, 0, [[0, 1, 2, 3]])
    assert int(dist_m4ri.DistanceResult(3, 6, 100)) == 6


def test_cache_newer_version_not_overwritten(tmp_path, capsys):
    import json
    cache_file = tmp_path / "newer_cache.json"
    newer = {"__version__": "99.0.0", "some_key": {"dist": 7}}
    cache_file.write_text(json.dumps(newer))
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    try:
        dist_m4ri.load_distance_cache(str(cache_file))
        assert "has newer version 99.0.0" in capsys.readouterr().err
        assert "some_key" not in dist_m4ri._distance_cache
        # neither a calculation with this cache file nor save_distance_cache() overwrites it
        d = dist_m4ri.compute_classical_distance(S5_H, method=2, wmax=3, threads=2, cache_file=str(cache_file))
        assert d == 2
        dist_m4ri.save_distance_cache(str(cache_file))
        assert json.loads(cache_file.read_text()) == newer
        # once the file is deleted with clear_distance_cache(clear_file=True), it is written again
        dist_m4ri.clear_distance_cache(str(cache_file), clear_file=True)
        dist_m4ri._distance_cache["code_test"] = {"dist": 3, "dmin": 3, "dmax": 3}
        dist_m4ri.save_distance_cache(str(cache_file))
        assert json.loads(cache_file.read_text())["__version__"] == dist_m4ri.__version__
    finally:
        dist_m4ri.clear_distance_cache()


def test_run_dist_m4ri_stop_event(capsys):
    import threading
    import time
    s5 = dict(finH=S5_H, finL=S5_L)
    # More than 64 KB of stderr (the matrices: debug=128|2048) with a stop_event which is not set: no deadlock
    result = {}

    def run():
        result["res"] = dist_m4ri.run_dist_m4ri(
            method=2, wmax=3, threads=2, debug=128 | 2048, stop_event=threading.Event(), **s5
        )

    worker = threading.Thread(target=run, daemon=True)
    worker.start()
    worker.join(timeout=120)
    assert not worker.is_alive(), "run_dist_m4ri() with stop_event is blocked by a full stderr pipe"
    assert tuple(result["res"]) == (4, 0, 0)
    assert len(capsys.readouterr().out) > 65536  # the stderr of the binary is printed (debug > 0)

    # Cancelled run: the binary is terminated, and RuntimeError is raised
    stop = threading.Event()
    timer = threading.Timer(0.3, stop.set)
    timer.start()
    t0 = time.time()
    try:
        with pytest.raises(RuntimeError, match="cancelled by stop_event"):
            dist_m4ri.run_dist_m4ri(method=1, steps=10**9, min_hits=0, timeout=60, threads=2, stop_event=stop, **s5)
        assert time.time() - t0 < 30
    finally:
        timer.cancel()


def test_add_noise_classically_controlled():
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circuit = stim.Circuit("""
        R 0 1 2
        TICK
        H 0
        CX 0 1
        TICK
        M 0
        CX rec[-1] 1
        CZ sweep[0] 2
        TICK
        M 1 2
        DETECTOR rec[-2]
        OBSERVABLE_INCLUDE(0) rec[-1]
    """)
    noisy = dist_m4ri.add_noise(circuit, p=0.01)
    assert "DEPOLARIZE2(0.01) 0 1" in str(noisy)
    # noise only on qubits: a classically controlled gate gets DEPOLARIZE1(p/10) on its qubit only
    for inst in noisy.flattened():
        if inst.name in ("DEPOLARIZE1", "DEPOLARIZE2"):
            assert all(t.qubit_value is not None for t in inst.targets_copy())
    lines = [s.strip() for s in str(noisy).splitlines()]
    assert lines[lines.index("CX rec[-1] 1") + 1] == "DEPOLARIZE1(0.001) 1"
    assert lines[lines.index("CZ sweep[0] 2") + 1] == "DEPOLARIZE1(0.001) 2"
    assert noisy.detector_error_model().num_errors > 0


def test_add_noise_gate_classes_and_strip_noise():
    # add_noise() classifies the gates with stim.gate_data(): DEPOLARIZE1(p/10) after every one-qubit gate,
    # DEPOLARIZE2(p) after every two-qubit gate, flips at resets and measurements, and no noise for the multi-qubit
    # measurements and Pauli product gates; strip_noise() removes all noise (heralded noise -> MPAD)
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circuit = stim.Circuit("""
        R 0 1 2 3
        RX 4
        TICK
        H_XY 0
        ISWAP 1 2
        C_XYZ 3
        SQRT_XX 0 4
        TICK
        MPP X1*X2
        MXX 1 2
        SPP Z3
        MRX 4
        MY 0
        TICK
        M 1 2 3
        DETECTOR rec[-2]
        OBSERVABLE_INCLUDE(0) rec[-1]
    """)
    noisy = dist_m4ri.add_noise(circuit, p=0.01)
    lines = [s.strip() for s in str(noisy).splitlines()]

    def after(line: str) -> str:
        return lines[lines.index(line) + 1]

    def before(line: str) -> str:
        return lines[lines.index(line) - 1]

    assert after("R 0 1 2 3") == "X_ERROR(0.01) 0 1 2 3" and after("RX 4") == "Z_ERROR(0.01) 4"
    assert after("H_XY 0") == "DEPOLARIZE1(0.001) 0" and after("C_XYZ 3") == "DEPOLARIZE1(0.001) 3"
    assert after("ISWAP 1 2") == "DEPOLARIZE2(0.01) 1 2" and after("SQRT_XX 0 4") == "DEPOLARIZE2(0.01) 0 4"
    assert after("MPP X1*X2") == "MXX 1 2" and after("MXX 1 2") == "SPP Z3"  # (no noise)
    assert before("MRX 4") == after("MRX 4") == "Z_ERROR(0.01) 4"
    assert before("MY 0") == "X_ERROR(0.01) 0" and before("M 1 2 3") == "X_ERROR(0.01) 1 2 3"
    assert not dist_m4ri.has_noise(circuit) and dist_m4ri.has_noise(noisy)
    # strip_noise() and force=True: the same circuit, and the same noise model
    assert dist_m4ri.strip_noise(noisy) == circuit
    assert dist_m4ri.add_noise(noisy, p=0.01, force=True) == noisy
    assert dist_m4ri.add_noise(noisy, p=0.0, force=True) == circuit and dist_m4ri.add_noise(circuit, p=0.0) == circuit
    # a noisy measurement is noise; heralded noise becomes MPAD (the measurement record does not change)
    other = stim.Circuit("R 0 1\nHERALDED_ERASE(0.01) 0\nM(0.02) 0 1\nDETECTOR rec[-1] rec[-3]\n")
    assert dist_m4ri.has_noise(stim.Circuit("R 0\nM(0.01) 0\n"))
    stripped = dist_m4ri.strip_noise(other)
    assert str(stripped).splitlines() == ["R 0 1", "MPAD 0", "M 0 1", "DETECTOR rec[-1] rec[-3]"]
    assert stripped.num_measurements == other.num_measurements and not dist_m4ri.has_noise(stripped)


def test_css_outc_sector_files_and_finc(tmp_path, capsys):
    cache_file = str(tmp_path / "css_cache.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    kw = dict(Hx=TRY_X, Hz=TRY_Z, method=2, wmax=8, threads=2)
    try:
        # A calculation writes the sector files outC_X / outC_Z ...
        res1 = dist_m4ri.compute_css_distance(outC=str(tmp_path / "cws1.nz"), cache_file=cache_file, **kw)
        assert res1 == (4, [4, 4, 0], [4, 4, 0])
        # ... and so does the cached exact result
        capsys.readouterr()
        res2 = dist_m4ri.compute_css_distance(
            outC=str(tmp_path / "cws2.nz"), cache_file=cache_file, debug=dist_m4ri.PY_DBG_CACHE, **kw
        )
        assert res2 == res1
        assert "Cache hit for CSS distance" in capsys.readouterr().out
        assert not (tmp_path / "cws2.nz").exists()
        for s in ("X", "Z"):
            cws1 = dist_m4ri.read_sparse_vectors(str(tmp_path / f"cws1_{s}.nz"))
            cws2 = dist_m4ri.read_sparse_vectors(str(tmp_path / f"cws2_{s}.nz"))
            assert len(cws1) == 40 and sorted(map(tuple, cws1)) == sorted(map(tuple, cws2))

        dist_m4ri.disable_distance_cache()
        # finC identical to outC: the sector files are read (finC itself does not exist)
        same = str(tmp_path / "cws1.nz")
        res3 = dist_m4ri.compute_css_distance(finC=same, outC=same, debug=1, **kw)
        assert res3 == res1
        out = capsys.readouterr().out
        assert "codewords from " + str(tmp_path / "cws1_X.nz") in out
        assert "codewords from " + str(tmp_path / "cws1_Z.nz") in out
        # Without sector files, finC itself is given to both runs (codewords not valid in a sector are skipped)
        mixed = tmp_path / "mixed.nz"
        mixed.write_text((tmp_path / "cws1_X.nz").read_text())
        res4 = dist_m4ri.compute_css_distance(finC=str(mixed), debug=1, **kw)
        assert res4 == res1
        assert capsys.readouterr().out.count("codewords from " + str(mixed)) == 2
    finally:
        dist_m4ri.enable_distance_cache()
        dist_m4ri.clear_distance_cache()


def _hgp_rep(a: int, b: int):
    """Hypergraph product of the repetition codes of lengths a and b (a planar surface code): dZ = a, dX = b."""
    def rep(n):
        H = np.zeros((n - 1, n), dtype=np.int8)
        for i in range(n - 1):
            H[i, i] = H[i, i + 1] = 1
        return H
    H1, H2 = rep(a), rep(b)
    Hx = np.hstack([np.kron(H1, np.eye(b, dtype=np.int8)), np.kron(np.eye(a - 1, dtype=np.int8), H2.T)]) % 2
    Hz = np.hstack([np.kron(np.eye(a, dtype=np.int8), H2), np.kron(H1.T, np.eye(b - 1, dtype=np.int8))]) % 2
    return Hx.astype(np.int8), Hz.astype(np.int8)


def test_dstop_stop_target(tmp_path, capsys):
    import time
    s5 = dict(finH=S5_H, finL=S5_L)  # d = 5
    # method=3: the run ends once CC has certified dmin >= dstop, without waiting for RW (here: until the timeout)
    t0 = time.time()
    res = dist_m4ri.run_dist_m4ri(method=3, dstop=4, steps=10**8, min_hits=0, timeout=60, threads=4, **s5)
    assert res[0] == 4 and res[1] in (0, 5) and time.time() - t0 < 30
    # method=2 needs no wmax with dstop; dstop is not an upper bound (dmax = 0: no codeword found)
    assert tuple(dist_m4ri.run_dist_m4ri(method=2, dstop=4, threads=2, **s5)) == (4, 0, 0)
    with pytest.raises(ValueError, match="dstop>0"):
        dist_m4ri.run_dist_m4ri(method=2, timeout=0, threads=2, **s5)

    cache_file = str(tmp_path / "dstop_cache.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    kw = dict(method=2, threads=2, return_info=True, cache_file=cache_file)
    try:
        assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, dstop=4, **kw) == (4, [4, 0, 0])
        # a cached lower bound dmin >= dstop answers the request
        capsys.readouterr()
        res2 = dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, dstop=3, debug=dist_m4ri.PY_DBG_CACHE, **kw)
        assert res2 == (4, [4, 0, 0])
        assert "Cache hit for quantum distance (dmin=4 >= dstop=3 known)!" in capsys.readouterr().out
        # a larger stop target continues from the cached lower bound, and so does the exact calculation
        assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, dstop=5, **kw) == (5, [5, 0, 0])
        assert dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, wmax=5, **kw) == (5, [5, 5, 0])
    finally:
        dist_m4ri.clear_distance_cache()

    text = dist_m4ri.explain_bounds([4, 0, 0], method=2, dstop=4)
    assert "Stop target (dstop = 4): reached" in text
    assert "Stop target" not in dist_m4ri.explain_bounds([4, 0, 0], method=2)
    assert "Stop target" not in dist_m4ri.explain_bounds([5, 5, 0], method=2, dstop=4)


def test_cli_dstop(capsys):
    assert dist_m4ri.parse_cli_args(["dstop=4"])["dstop"] == 4
    assert dist_m4ri.parse_cli_args([])["dstop"] == 0
    try:
        assert dist_m4ri.main([f"finH={S5_H}", "dstop=-1", "--no-cache"]) == 1
        assert "invalid 'dstop=-1': dstop must be non-negative" in capsys.readouterr().err
        ret = dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "method=2", "dstop=4", "threads=2", "--no-cache", "-v"])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.rstrip().endswith("4 0 0") and "Stop target (dstop = 4): reached" in out
        # CSS: a user dmax (a bound on d = min(dX, dZ)) only stops both sector runs; then d = dmax
        ret = dist_m4ri.main([f"Hx={TRY_X}", f"Hz={TRY_Z}", "method=2", "dmax=4", "threads=2", "--no-cache"])
        assert ret == 0
        assert "dX: 4 0 0  dZ: 4 0 0  (d = 4) (exact)" in capsys.readouterr().out
    finally:
        dist_m4ri.enable_distance_cache()


def test_css_user_dmax_and_sector_records(tmp_path, capsys):
    Hx, Hz = _hgp_rep(3, 5)  # dZ = 3, dX = 5
    cache_file = str(tmp_path / "css_sectors.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    kw = dict(method=2, wmax=6, threads=2, cache_file=cache_file)
    try:
        # A user dmax is a bound on d = min(dX, dZ), not on either sector: it stops both sector runs (dstop),
        # and it is neither reported as nor stored as the upper bound of a sector
        assert dist_m4ri.compute_css_distance(Hx, Hz, dmax=3, **kw) == (3, [3, 0, 0], [3, 0, 0])
        assert dist_m4ri.get_cached_distance(H=Hz, G=Hx, cache_file=cache_file)["d_info"] == [3, 0, 0]
        # Without dmax, both sectors are computed (from the cached lower bounds)
        assert dist_m4ri.compute_css_distance(Hx, Hz, **kw) == (3, [5, 5, 0], [3, 3, 0])
        # The sector records are those of compute_quantum_distance(): (H=Hz, G=Hx) -> dX, (H=Hx, G=Hz) -> dZ
        capsys.readouterr()
        dbg = dist_m4ri.PY_DBG_CACHE
        assert dist_m4ri.compute_quantum_distance(Hz, G=Hx, return_info=True, debug=dbg, **kw) == (5, [5, 5, 0])
        assert dist_m4ri.compute_quantum_distance(Hx, G=Hz, return_info=True, debug=dbg, **kw) == (3, [3, 3, 0])
        assert capsys.readouterr().out.count("Cache hit for quantum distance (exact distance known)!") == 2
        entry = dist_m4ri.get_cached_distance(Hx=Hx, Hz=Hz, cache_file=cache_file)
        assert entry["dX"] == [5, 5, 0] and entry["dZ"] == [3, 3, 0]
        assert (entry["dist"], entry["dmin"], entry["dmax"]) == (3, 3, 3)
        assert dist_m4ri.compute_css_distance(Hx, Hz, dstop=3, debug=dbg, **kw) == (3, [5, 5, 0], [3, 3, 0])
        assert "Cache hit for CSS distance (exact distance known)!" in capsys.readouterr().out
    finally:
        dist_m4ri.clear_distance_cache()


def test_css_cache_migration(tmp_path):
    import json
    Hx, Hz = _hgp_rep(3, 5)  # dZ = 3, dX = 5
    st = dist_m4ri.get_sparse_array_state
    key = f"css:X={st(Hx)}:Z={st(Hz)}"
    key_x, key_z = f"quantum:H={st(Hz)}:G={st(Hx)}", f"quantum:H={st(Hx)}:G={st(Hz)}"
    # A CSS record of an earlier version after compute_css_distance(dmax=3): the bound on d = min(dX, dZ) was
    # reported back as the upper bound of both sectors (wrong for dX = 5)
    old = {"dist": 3, "dmin": 3, "dmax": 3, "rw_steps": 0, "dmin_X": 3, "dmax_X": 3, "rw_steps_X": 0,
           "dmin_Z": 3, "dmax_Z": 3, "rw_steps_Z": 0, "dX": [3, 3, 0], "dZ": [3, 3, 0], "cws_X": [], "cws_Z": []}
    cache_file = tmp_path / "old_css_cache.json"
    cache_file.write_text(json.dumps({"__version__": dist_m4ri.__version__, key: old}))
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    try:
        # The lower bounds are kept, the unconfirmed upper bounds are dropped, and the old record is removed
        entry = dist_m4ri.get_cached_distance(Hx=Hx, Hz=Hz, cache_file=str(cache_file))
        assert entry["dX"] == [3, 0, 0] and entry["dZ"] == [3, 0, 0]
        data = json.loads(cache_file.read_text())
        assert key not in data and key_x in data and key_z in data
        res = dist_m4ri.compute_css_distance(Hx, Hz, method=2, wmax=6, threads=2, cache_file=str(cache_file))
        assert res == (3, [5, 5, 0], [3, 3, 0])

        # An upper bound confirmed by a stored codeword (only its weight is used here), or smaller than the upper
        # bound of the other sector, is kept
        dist_m4ri.clear_distance_cache()
        old2 = dict(old, dmax_X=5, cws_X=[[0, 1, 2, 3, 4]], dmax_Z=3)
        cache_file.write_text(json.dumps({"__version__": dist_m4ri.__version__, key: old2}))
        entry2 = dist_m4ri.get_cached_distance(Hx=Hx, Hz=Hz, cache_file=str(cache_file))
        assert entry2["dX"] == [3, 5, 0] and entry2["dZ"] == [3, 3, 0] and entry2["cws_X"] == [[0, 1, 2, 3, 4]]
        assert (entry2["dist"], entry2["dmin"], entry2["dmax"]) == (3, 3, 3)
    finally:
        dist_m4ri.clear_distance_cache()


def test_cached_outc_keeps_finc_codewords(tmp_path):
    cache_file = str(tmp_path / "cws_cache.json")
    out = str(tmp_path / "cws.nz")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    kw = dict(method=2, wmax=5, threads=2, cache_file=cache_file)
    try:
        d, cws = dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, outC=out, do_cws=True, **kw)
        assert d == 5 and len(cws) >= 2
        extra = sorted(set(cws[0]) ^ set(cws[1]))  # another vector (of even weight, not in cws)
        dist_m4ri._write_nzlist_file(out, cws + [extra])
        # A cache hit with finC identical to outC keeps the codewords in the file
        d2, cws2 = dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, finC=out, outC=out, do_cws=True, **kw)
        assert d2 == 5 and sorted(map(tuple, cws2)) == sorted(map(tuple, cws))
        assert sorted(map(tuple, dist_m4ri.read_sparse_vectors(out))) == sorted(map(tuple, cws + [extra]))
        # Without finC, outC gets the cached codewords
        dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, outC=out, **kw)
        assert sorted(map(tuple, dist_m4ri.read_sparse_vectors(out))) == sorted(map(tuple, cws))
    finally:
        dist_m4ri.clear_distance_cache()


SURF_D3_DEM = os.path.join(EXAMPLES_DIR, "surf_d3.dem")


def test_unknown_keyword_arguments():
    with pytest.raises(TypeError, match=r"compute_quantum_distance\(\) got an unexpected keyword argument 'dmx' "
                                        r"\(did you mean 'dmax'\?\)"):
        dist_m4ri.compute_quantum_distance(S5_H, L=S5_L, method=2, wmax=5, dmx=5)
    with pytest.raises(TypeError, match=r"compute_css_distance\(\) got unexpected keyword arguments 'stpes'.*, "
                                        r"'dexp' \(did you mean 'd_exp'\?\)"):
        dist_m4ri.compute_css_distance(TRY_X, TRY_Z, method=2, wmax=4, stpes=10, dexp=4)
    for func, kw in ((dist_m4ri.compute_classical_distance, dict(H=S5_H)),
                     (dist_m4ri.compute_dem_distance, dict(dem=SURF_D3_DEM))):
        with pytest.raises(TypeError, match="unexpected keyword argument 'foo'"):
            func(foo=1, method=2, wmax=3, **kw)
    # the alias win (of kwin) is accepted
    assert dist_m4ri.compute_classical_distance(S5_H, method=2, wmax=3, threads=2, win=0) == 2


def test_cli_conflict_warnings(capsys):
    try:
        ret = dist_m4ri.main([f"fdem={SURF_D3_DEM}", f"finH={S5_H}", "--simple", "method=2", "wmax=3", "threads=2",
                              "--no-cache"])
        assert ret == 0
        captured = capsys.readouterr()
        assert f"# Warning: fdem={SURF_D3_DEM} given: ignoring finH=" in captured.err
        assert "# Warning: --simple: disabled since version 0.13.0, ignored" in captured.err
        assert captured.out.rstrip().endswith("3 3 0 (exact)")

        ret = dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", f"Lx={S5_L}", "pmin=0.01", "classical=1", "method=2",
                              "wmax=5", "threads=2", "--no-cache"])
        assert ret == 0
        err = capsys.readouterr().err
        assert "# Warning: pmin=: only used with fdem= (a DEM or Stim circuit); ignored" in err
        assert "# Warning: Lx=: only used for a CSS code (Hx= and Hz=); ignored" in err
        assert "# Warning: classical=1 ignored: the input is a quantum code or a DEM" in err

        prefix = os.path.join(EXAMPLES_DIR, "try")
        ret = dist_m4ri.main([f"fin={prefix}", f"finH={TRY_X}", "method=2", "wmax=4", "threads=2", "--no-cache"])
        assert ret == 0
        captured = capsys.readouterr()
        assert f"# Warning: fin={prefix} together with finH=: the explicit file names take precedence" in captured.err
        assert captured.out.rstrip().endswith("4 4 0 (exact)")

        ret = dist_m4ri.main([f"Hx={TRY_X}", f"Hz={TRY_Z}", f"finH={S5_H}", "method=2", "wmax=4", "threads=2",
                              "--no-cache"])
        assert ret == 0
        assert "# Warning: Hx=/Hz= given (CSS code): ignoring finH=" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()


def test_cli_outc_keeps_finc_codewords(tmp_path, capsys):
    cache_file = str(tmp_path / "cli_cache.json")
    out = str(tmp_path / "cli_cws.nz")
    base = [f"finH={S5_H}", f"finL={S5_L}", "method=2", "wmax=5", "threads=2", f"cache={cache_file}"]
    dist_m4ri.clear_distance_cache()
    try:
        assert dist_m4ri.main(base + [f"outC={out}"]) == 0
        cws = dist_m4ri.read_sparse_vectors(out)
        assert len(cws) >= 2
        extra = sorted(set(cws[0]) ^ set(cws[1]))
        dist_m4ri._write_nzlist_file(out, cws + [extra])
        # a cache hit with finC identical to outC: the file is written once, keeping its codewords
        assert dist_m4ri.main(base + [f"finC={out}", f"outC={out}"]) == 0
        assert sorted(map(tuple, dist_m4ri.read_sparse_vectors(out))) == sorted(map(tuple, cws + [extra]))
        capsys.readouterr()
    finally:
        dist_m4ri.clear_distance_cache()


def test_dem_cache_merged_bounds_and_circuit_key(tmp_path, capsys):
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circ = stim.Circuit.generated(
        "surface_code:rotated_memory_z", rounds=3, distance=3, after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001, after_reset_flip_probability=0.001,
    )
    cache_file = str(tmp_path / "dem_cache.json")
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    kw = dict(method=2, threads=2, cache_file=cache_file)
    try:
        assert dist_m4ri.get_cached_distance(circuit=circ, cache_file=cache_file) is None
        assert dist_m4ri.compute_dem_distance(circuit=circ, dstop=2, **kw) == (2, [2, 0, 0])
        # get_cached_distance(circuit=...) finds the record of the processed DEM (with all detectors)
        assert dist_m4ri.get_cached_distance(circuit=circ, cache_file=cache_file)["d_info"] == [2, 0, 0]
        assert dist_m4ri.get_cached_distance(circuit=circ, rounds=4, cache_file=cache_file) is None
        # full= is disabled (ignored with a warning): the same record
        assert dist_m4ri.get_cached_distance(circuit=circ, full=True, cache_file=cache_file)["d_info"] == [2, 0, 0]
        assert "# Warning: full=: disabled since version 0.13.0" in capsys.readouterr().err
        # the exact distance from the cached lower bound; the cached dist is that of the merged bounds
        d, info, cws = dist_m4ri.compute_dem_distance(circuit=circ, wmax=3, do_cws=True, **kw)
        assert d == 3 and info == [3, 3, 0] and len(cws) > 0
        entry = dist_m4ri.get_cached_distance(circuit=circ, cache_file=cache_file)
        assert entry["dist"] == 3 and entry["d_info"] == [3, 3, 0] and len(entry["cws"]) == len(cws)
        stim_file = tmp_path / "c.stim"
        circ.to_file(str(stim_file))
        assert dist_m4ri.get_cached_distance(circuit=str(stim_file), cache_file=cache_file)["dist"] == 3
    finally:
        dist_m4ri.clear_distance_cache()


def test_dem_merge_equal_columns(tmp_path, capsys):
    # the binary merges the equal columns of a DEM (merge=True, the default; exact): the codewords use the first one of
    # each set of equal columns (0-based columns 0=4, 1=3, and 5=6); merge=False (CLI: merge=0) keeps all columns
    dem_file = str(tmp_path / "dup.dem")
    with open(dem_file, "w") as f:
        f.write("error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\nerror(0.1) D0 D1\nerror(0.1) D0 L0\n"
                "error(0.1) D2\nerror(0.1) D2\n")
    dist_m4ri.disable_distance_cache()
    try:
        kw = dict(dem=dem_file, method=2, wmax=3, do_cws=True, threads=1)
        assert dist_m4ri.compute_dem_distance(**kw) == (3, [3, 3, 0], [[0, 1, 2]])
        d, info, cws = dist_m4ri.compute_dem_distance(merge=False, **kw)
        assert (d, info) == (3, [3, 3, 0])
        assert sorted(map(tuple, cws)) == [(0, 1, 2), (0, 2, 3), (1, 2, 4), (2, 3, 4)]
        out = str(tmp_path / "cws.nz")
        assert dist_m4ri.main([f"fdem={dem_file}", "merge=0", "method=2", "wmax=3", "threads=1", f"outC={out}",
                               "--no-cache"]) == 0
        assert len(dist_m4ri.read_sparse_vectors(out)) == 4
        capsys.readouterr()
        assert dist_m4ri.main([f"fdem={dem_file}", "merge=2", "method=2", "wmax=3", "--no-cache"]) == 1
        assert "invalid 'merge=2': merge should be 0 or 1" in capsys.readouterr().err
        assert dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", "merge=0", "method=2", "wmax=5", "threads=2",
                               "--no-cache"]) == 0
        assert "# Warning: merge=: only used with fdem= (a DEM or Stim circuit); ignored" in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()


def test_dem_expand_equal_columns(tmp_path, capsys):
    # expand=True (Python only): each codeword of the merged list is listed for every choice among equal columns
    # (0-based columns 0=4, 1=3, and 5=6); the cache keeps the merged list; maxC limits the expanded list (warning)
    dem_file = str(tmp_path / "dup.dem")
    with open(dem_file, "w") as f:
        f.write("error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\nerror(0.1) D0 D1\nerror(0.1) D0 L0\n"
                "error(0.1) D2\nerror(0.1) D2\n")
    full = [[0, 1, 2], [0, 2, 3], [1, 2, 4], [2, 3, 4]]
    cache_file = str(tmp_path / "cache.json")
    kw = dict(dem=dem_file, method=2, wmax=3, do_cws=True, threads=1, cache_file=cache_file)
    dist_m4ri.clear_distance_cache()
    dist_m4ri.enable_distance_cache()
    try:
        assert dist_m4ri.compute_dem_distance(**kw) == (3, [3, 3, 0], [[0, 1, 2]])
        # from the cache (the merged list), expanded
        d, info, cws = dist_m4ri.compute_dem_distance(expand=True, **kw)
        assert (d, info) == (3, [3, 3, 0]) and sorted(cws) == full
        assert dist_m4ri.get_cached_distance(dem=dem_file, cache_file=cache_file)["cws"] == [[0, 1, 2]]
        capsys.readouterr()
        d, info, cws = dist_m4ri.compute_dem_distance(expand=True, maxC=2, **kw)
        assert len(cws) == 2 and all(cw in full for cw in cws)
        assert ("# Warning: expand: the expanded list has 4 codewords (1 up to equal columns); only the first maxC=2 "
                "are written (the list is incomplete)") in capsys.readouterr().err
        # a run of the binary (no cache), CLI
        out = str(tmp_path / "cws.nz")
        assert dist_m4ri.main([f"fdem={dem_file}", "expand=1", "method=2", "wmax=3", "threads=1", f"outC={out}",
                               "--no-cache"]) == 0
        assert sorted(dist_m4ri.read_sparse_vectors(out)) == full
        capsys.readouterr()
        assert dist_m4ri.main([f"fdem={dem_file}", "expand=2", "method=2", "wmax=3", "--no-cache"]) == 1
        assert "invalid 'expand=2': expand should be 0 or 1" in capsys.readouterr().err
        assert dist_m4ri.main([f"fdem={dem_file}", "expand=1", "method=2", "wmax=3", "threads=1", "--no-cache"]) == 0
        assert "# Warning: expand=: only used with outC= or --cws; ignored" in capsys.readouterr().err
    finally:
        dist_m4ri.clear_distance_cache()
        dist_m4ri.enable_distance_cache()


def test_dem_error_columns_parser(tmp_path):
    # the columns of a DEM as read by the binary: repeat blocks unrolled with shift_detectors, components summed mod 2,
    # and the errors with probability < pmin skipped
    dem_file = tmp_path / "r.dem"
    dem_file.write_text("error(0.1) D0 L0\nrepeat 2 {\n    error(0.1) D0 D1\n    shift_detectors 1\n}\n"
                        "error(0.01) D5\nerror(0.1) D0 ^ D0 D2  # comment\ndetector(1, 2) D3\n")
    cols = dist_m4ri._dem_error_columns(str(dem_file))
    assert cols == [({0}, {0}), ({0, 1}, set()), ({1, 2}, set()), ({7}, set()), ({4}, set())]
    assert dist_m4ri._dem_error_columns(str(dem_file), pmin=0.05) == [cols[i] for i in (0, 1, 2, 4)]
    if dist_m4ri._HAS_STIM:  # the same columns from the stim object (flattened)
        import stim
        assert dist_m4ri._dem_error_columns(stim.DetectorErrorModel(dem_file.read_text())) == cols


def test_dem_circuit_outc_without_out_dem_warning(tmp_path, capsys):
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circ = stim.Circuit.generated("repetition_code:memory", rounds=2, distance=3, after_clifford_depolarization=0.01)
    kw = dict(circuit=circ, method=2, wmax=3, threads=1, outC=str(tmp_path / "cws.nz"), out_dir=str(tmp_path))
    dist_m4ri.disable_distance_cache()
    try:
        assert dist_m4ri.compute_dem_distance(**kw)[0] == 3
        assert "the codewords index the error lines of the DEM generated from the Stim circuit in this run" in \
            capsys.readouterr().err
        assert dist_m4ri.compute_dem_distance(out_dem=str(tmp_path / "c.dem"), **kw)[0] == 3
        assert "the codewords index the error lines" not in capsys.readouterr().err
    finally:
        dist_m4ri.enable_distance_cache()


def test_dem_drop_minority_detectors(tmp_path, capsys):
    # drop=auto (the default): the binary finds the detectors which can be dropped exactly (all errors touching them
    # have twins), here the minority-basis detectors of a CSS circuit, with the distance and codewords of the full DEM;
    # drop=css: the minority-basis detectors of the classification of the circuit, checked by the binary
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circ = stim.Circuit.generated(
        "surface_code:rotated_memory_z", rounds=3, distance=3, after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001, after_reset_flip_probability=0.001,
    )
    drop = []
    dem, _, _ = dist_m4ri._dem_from_input(None, circ, drop_out=drop)
    assert len(drop) == 8 and dem.num_detectors == 24  # the X detectors of rounds 2 and 3 (4 each)
    dist_m4ri.disable_distance_cache()
    try:
        out_d = tmp_path / "d.txt"
        kw = dict(circuit=circ, method=2, wmax=3, threads=1, do_cws=True, outD=str(out_d))
        cols = dist_m4ri._dem_error_columns(dem)
        for mode in (None, "auto", 1, True, "css"):
            d, info, cws = dist_m4ri.compute_dem_distance(debug=1, drop=mode, **kw)
            out = capsys.readouterr().out
            assert (d, info) == (3, [3, 3, 0]) and len(cws) > 0
            assert re.search(r"dropped 8 of 24 detectors \(drop=.*\), exact: the \d+ columns touching them have twins "
                             r"\(the same observables\)", out)
            # the detectors dropped (outD): those of the classification, also with drop=auto
            assert dist_m4ri.read_drop_file(out_d) == (drop, 24)
            label = "css" if mode == "css" else "auto"
            assert out_d.read_text().startswith(f"# 8 of 24 detectors dropped (drop={label}), 0-based\n")
            for cw in cws:  # codewords of the full DEM
                dets, obs = set(), set()
                for j in cw:
                    dets ^= cols[j][0]
                    obs ^= cols[j][1]
                assert not dets and obs
        for mode in (False, 0, "0"):
            assert dist_m4ri.compute_dem_distance(drop=mode, debug=1, **kw)[:2] == (3, [3, 3, 0])
            assert "dropped" not in capsys.readouterr().out
            assert dist_m4ri.read_drop_file(out_d) == ([], 24)
        with pytest.raises(ValueError, match="no such file"):
            dist_m4ri.compute_dem_distance(circuit=circ, drop=str(tmp_path / "none.txt"))
        # a DEM file with the detectors to drop in a file: A != 0 (exact), and a column without a twin (a note)
        dem_file, drop_file = tmp_path / "a.dem", tmp_path / "drop.txt"
        drop_file.write_text("# the other basis\n2\n")
        dem_file.write_text("error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\nerror(0.1) D0 D2\nerror(0.1) D2 L0\n"
                            "error(0.1) D1 D2 L0\n")
        kw = dict(dem=str(dem_file), drop=str(drop_file), method=2, wmax=3, threads=1)
        assert dist_m4ri.compute_dem_distance(**kw) == (3, [3, 3, 0])
        assert "# Note: dropping" not in capsys.readouterr().err
        with open(dem_file, "a") as f:
            f.write("error(0.1) D2 D3\n")
        assert dist_m4ri.compute_dem_distance(**kw) == (3, [3, 3, 0])
        assert (f"# Note: dropping the detectors in {drop_file} is not exact for this DEM (column 7 (1-based) "
                f"touches a dropped detector") in capsys.readouterr().err
        drop_file.write_text("9\n")  # (not a detector of the DEM)
        assert dist_m4ri.compute_dem_distance(**kw) == (3, [3, 3, 0])
        assert f"# Note: none of the detectors in {drop_file} is in the DEM; nothing dropped" in capsys.readouterr().err
        # CLI: the values of drop, drop=css for a DEM file (nothing dropped), a missing drop file (an error), drop=,
        # outD=, and --out-reduced without fdem, and basis= without drop=css (warnings)
        for val, exp in (("1", True), ("auto", True), ("true", True), ("0", False), ("off", False), ("css", "css"),
                         ("CSS", "css"), ("d.txt", "d.txt")):
            assert dist_m4ri.parse_cli_args([f"drop={val}"])["drop"] == exp
        assert dist_m4ri.main([f"fdem={dem_file}", "drop=css", "method=2", "wmax=3", "threads=1", "--no-cache"]) == 0
        assert "# Warning: drop=css: no minority-basis detectors" in capsys.readouterr().err
        assert dist_m4ri.compute_dem_distance(dem=str(dem_file), drop="css", method=2, wmax=3, threads=1,
                                              outD=str(out_d)) == (3, [3, 3, 0])
        assert "using drop=auto" in capsys.readouterr().err
        assert out_d.read_text().startswith("# ") and "(drop=auto" in out_d.read_text()
        assert dist_m4ri.main([f"fdem={dem_file}", f"drop={tmp_path / 'none.txt'}", "method=2", "wmax=3",
                               "--no-cache"]) == 1
        assert "none.txt: no such file" in capsys.readouterr().err
        assert dist_m4ri.main([f"finH={S5_H}", f"finL={S5_L}", f"drop={drop_file}", f"outD={out_d}", "--out-reduced",
                               "method=2", "wmax=5", "threads=2", "--no-cache"]) == 0
        assert ("# Warning: drop=, outD=, --out-reduced: only used with fdem= (a DEM or Stim circuit); ignored" in
                capsys.readouterr().err)
        stim_file = tmp_path / "surf_d3_Z.stim"
        circ.to_file(str(stim_file))
        for opts, warn in (([], True), (["drop=css"], False), (["-v"], False)):
            assert dist_m4ri.main([str(stim_file), "basis=Z", "method=2", "wmax=3", "threads=1", "--no-cache"] +
                                  opts) == 0
            assert ("# Warning: basis=: only used with drop=css" in capsys.readouterr().err) == warn
    finally:
        dist_m4ri.enable_distance_cache()


def test_dem_out_reduced(tmp_path, capsys):
    # outD and out_reduced: the detectors dropped, and the reduced DEM with its maps to the DEM (also with a cached
    # result: the binary runs only for the detectors dropped)
    dem_file = tmp_path / "surf_d3_full.dem"
    shutil.copyfile(SURF_D3_DEM, dem_file)
    red_dem, red_json = tmp_path / "surf_d3_red.dem", tmp_path / "surf_d3_red.json"
    out_d = tmp_path / "d.txt"
    dist_m4ri.clear_distance_cache()
    try:
        for run in range(2):
            for f in (red_dem, red_json, out_d):
                if f.exists():
                    f.unlink()
            res = dist_m4ri.compute_dem_distance(dem=str(dem_file), method=2, wmax=3, threads=1, out_reduced=True,
                                                 outD=str(out_d), cache_file=str(tmp_path / "cache.json"), verbose=True)
            assert res == (3, [3, 3, 0])
            out = capsys.readouterr().out
            assert ("Cached result: running the binary only for the detectors dropped (outD)" in out) == (run == 1)
            assert f"Saved the reduced DEM (8 detectors dropped) to '{red_dem}'" in out
            assert red_dem.exists() and red_json.exists()
            assert dist_m4ri.read_drop_file(out_d) == ([5, 7, 8, 10, 13, 15, 16, 18], 24)
        info = json.loads(red_json.read_text())
        assert (info["__version__"], info["source"], info["dem"], info["drop"]) == \
            (dist_m4ri.__version__, str(dem_file), str(dem_file), "auto")
        assert (info["num_detectors"], info["num_errors"], info["pmin"]) == (24, 221, 0.0)
        assert info["dropped"] == [5, 7, 8, 10, 13, 15, 16, 18]
        assert info["detectors"] == [r for r in range(24) if r not in info["dropped"]]
        cols = dist_m4ri._dem_error_columns(str(dem_file))
        red_cols = dist_m4ri._dem_error_columns(str(red_dem))
        assert len(red_cols) == len(info["errors"]) == 55
        for (S, O), j in zip(red_cols, info["errors"]):  # the maps of the detectors and errors
            assert {info["detectors"][r] for r in S} == cols[j][0] and O == cols[j][1]
        assert all(not (S & set(info["dropped"])) for S, _ in (cols[j] for j in info["errors"]))
        # the reduced DEM has the same distance, and its codewords map to codewords of the DEM
        dist_m4ri.disable_distance_cache()
        d, d_info, cws = dist_m4ri.compute_dem_distance(dem=str(red_dem), method=2, wmax=3, threads=1, drop=0,
                                                        do_cws=True)
        assert (d, d_info, len(cws)) == (3, [3, 3, 0], 21)
        for cw in cws:
            dets, obs = set(), set()
            for i in cw:
                dets ^= cols[info["errors"][i]][0]
                obs ^= cols[info["errors"][i]][1]
            assert not dets and obs
        # CLI: --out-reduced PREFIX (a path), and write_reduced_dem() keeps the number of detectors and observables
        assert dist_m4ri.main([f"fdem={dem_file}", "--out-reduced", str(tmp_path / "x" / "pfx"), "method=2", "wmax=3",
                               "threads=1", "--no-cache"]) == 0
        assert (tmp_path / "x" / "pfx_red.dem").read_text().splitlines()[1:] == red_dem.read_text().splitlines()[1:]
        args = dist_m4ri.parse_cli_args(["--out-reduced", str(dem_file)])  # (an existing input file, not a prefix)
        assert (args["out_reduced"], args["fdem"]) == (True, str(dem_file))
        assert dist_m4ri.parse_cli_args(["out_reduced=0", "outD=y.txt"])["out_reduced"] is None
        small = tmp_path / "s.dem"
        small.write_text("error(0.1) D0 L0\nerror(0.2) D0 D1\nerror(0.3) D2 L1\n")
        dem_s, json_s = dist_m4ri.write_reduced_dem(str(small), [1, 2], str(tmp_path / "s"), num_detectors=4)
        assert Path(dem_s).read_text().splitlines()[1:] == ["error(0.1) D0 L0", "detector D1", "logical_observable L1"]
        assert json.loads(Path(json_s).read_text())["detectors"] == [0, 3]
    finally:
        dist_m4ri.clear_distance_cache()
        dist_m4ri.enable_distance_cache()


def test_dem_force_noise_and_p_noise(tmp_path, capsys):
    # p_noise: the strength of the noise added to a noiseless circuit (0: none); force_noise: the noise of the circuit
    # is replaced (p_noise=0: removed), e.g., the same noise model for circuits with other noise
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    c0 = stim.Circuit.generated("surface_code:rotated_memory_z", rounds=3, distance=3)
    cn = stim.Circuit.generated("surface_code:rotated_memory_z", rounds=3, distance=3,
                                after_clifford_depolarization=0.002, before_measure_flip_probability=0.003)
    noisy_file, cache = tmp_path / "surf_Z.stim", str(tmp_path / "cache.json")
    cn.to_file(str(noisy_file))
    kw = dict(method=2, wmax=3, threads=1)
    dist_m4ri.clear_distance_cache()
    try:
        # the same DEM (and cache record) as for the noiseless circuit with the noise added
        assert dist_m4ri.compute_dem_distance(circuit=cn, force_noise=True, out_dem=str(tmp_path / "f.dem"),
                                              cache_file=cache, **kw) == (3, [3, 3, 0])
        assert dist_m4ri.compute_dem_distance(circuit=c0, out_dem=str(tmp_path / "a.dem"), **kw) == (3, [3, 3, 0])
        assert (tmp_path / "f.dem").read_text() == (tmp_path / "a.dem").read_text()
        assert dist_m4ri.get_cached_distance(circuit=c0, cache_file=cache)["dmin"] == 3
        assert dist_m4ri.get_cached_distance(circuit=cn, force_noise=True, cache_file=cache)["dmin"] == 3
        assert dist_m4ri.get_cached_distance(circuit=cn, cache_file=cache) is None  # (its own noise: another DEM)
        # p_noise=0: no noise added (no distance for a noiseless circuit); with force_noise, all noise is removed
        dist_m4ri.disable_distance_cache()
        assert dist_m4ri.compute_dem_distance(circuit=c0, p_noise=0, **kw) == (0, [0, 0, 0])
        assert "the DEM has no error mechanisms" in capsys.readouterr().err
        stripped = tmp_path / "s.stim"
        assert dist_m4ri.compute_dem_distance(circuit=cn, p_noise=0, force_noise=True, out_stim=str(stripped),
                                              **kw) == (0, [0, 0, 0])
        assert stim.Circuit.from_file(str(stripped)) == c0
        assert dist_m4ri.compute_dem_distance(circuit=cn, p_noise=0, **kw) == (3, [3, 3, 0])  # (its own noise)
        with pytest.raises(ValueError, match="noise strength"):
            dist_m4ri.compute_dem_distance(circuit=c0, p_noise=0.7)
        # CLI: --p-noise P (p_noise=P), --force-noise (--force, force_noise=1); for a DEM file: ignored (a warning)
        for argv, exp in ((["--p-noise", "0.002"], (0.002, False)), (["p_noise=0"], (0.0, False)),
                          (["--force-noise"], (None, True)), (["--force"], (None, True)),
                          (["force_noise=1", "pnoise=0.01"], (0.01, True))):
            args = dist_m4ri.parse_cli_args(argv + [str(noisy_file)])
            assert (args["p_noise"], args["force_noise"], args["fdem"]) == exp + (str(noisy_file),)
        assert dist_m4ri.main([str(noisy_file), "--force-noise", "--p-noise", "0", "--out-stim", str(stripped),
                               "--no-cache"]) == 0
        out = capsys.readouterr()
        assert out.out.strip() == "0 0 0" and "no error mechanisms" in out.err
        assert dist_m4ri.main([str(noisy_file), "--force", "method=2", "wmax=3", "threads=1", "--no-cache"]) == 0
        assert capsys.readouterr().out.strip() == "3 3 0 (exact)"
        assert dist_m4ri.main([str(noisy_file), "p_noise=0.7", "--no-cache"]) == 1
        assert "the noise strength must be in [0, 0.5]" in capsys.readouterr().err
        assert dist_m4ri.main([f"fdem={SURF_D3_DEM}", "--force-noise", "p_noise=0.01", "method=2", "wmax=3",
                               "threads=1", "--no-cache"]) == 0
        assert ("# Warning: --p-noise, --force-noise: only used for a Stim circuit (.stim); ignored" in
                capsys.readouterr().err)
    finally:
        dist_m4ri.clear_distance_cache()
        dist_m4ri.enable_distance_cache()


def test_stim_pauli_channel_circuit(tmp_path, capsys):
    # an error channel with disjoint cases which stim cannot write as independent errors (here: PAULI_CHANNEL_1 with X
    # and Z, but no Y) needs approximate_disjoint_errors=True in stim
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    import stim
    circ = stim.Circuit.generated("surface_code:rotated_memory_z", rounds=3, distance=3,
                                  before_round_data_depolarization=0.001, before_measure_flip_probability=0.001)
    circ = stim.Circuit(str(circ).replace("DEPOLARIZE1(0.001)", "PAULI_CHANNEL_1(0.0005, 0, 0.0005)"))
    assert "PAULI_CHANNEL_1" in str(circ) and dist_m4ri.has_noise(circ)
    with pytest.raises(ValueError):
        circ.detector_error_model(decompose_errors=False)  # (without the approximation)
    dem, decomp, approx = dist_m4ri._circuit_dem(circ)
    assert approx is True and dem.num_errors > 0
    dist_m4ri.disable_distance_cache()
    try:
        assert dist_m4ri.compute_dem_distance(circuit=circ, method=2, wmax=3, threads=2) == (3, [3, 3, 0])
        stim_file = tmp_path / "pc.stim"
        circ.to_file(str(stim_file))
        assert dist_m4ri.main([str(stim_file), "method=2", "wmax=3", "threads=2", "--no-cache", "-v"]) == 0
        out = capsys.readouterr().out
        assert "approximate_disjoint_errors=True" in out and "3 3 0 (exact)" in out
        # a circuit without such channels is converted as before
        assert dist_m4ri._circuit_dem(stim.Circuit.generated("surface_code:rotated_memory_z", rounds=2, distance=3,
                                                             after_clifford_depolarization=0.001))[2] is False
    finally:
        dist_m4ri.enable_distance_cache()


def _surface_memory_z(x_in_h_frame: bool, meas_before_h: bool, rounds: int = 3):
    """
    Memory-Z circuit of the d=3 rotated surface code (data 0..8, Z ancillas 9..12, X ancillas 13..16).  With
    x_in_h_frame, the X checks are measured as Z checks in the Hadamard frame of the data qubits (H on all data, CX
    data -> ancilla, H on all data), with the X ancillas measured before or after the closing H gates.
    """
    import stim
    x_checks = [[1, 2, 4, 5], [3, 4, 6, 7], [0, 1], [7, 8]]
    z_checks = [[0, 1, 3, 4], [4, 5, 7, 8], [3, 6], [2, 5]]
    data, za, xa = list(range(9)), [9, 10, 11, 12], [13, 14, 15, 16]

    def line(name, qubits):
        return f"{name} " + " ".join(map(str, qubits)) + "\n"

    def body():
        s = "".join(line("CX", [x for q in row for x in (q, a)]) for a, row in zip(za, z_checks))
        if x_in_h_frame:
            s += line("MR", za) + line("H", data)
            s += "".join(line("CX", [x for q in row for x in (q, a)]) for a, row in zip(xa, x_checks))
            s += (line("MR", xa) + line("H", data)) if meas_before_h else (line("H", data) + line("MR", xa))
        else:
            s += line("H", xa)
            s += "".join(line("CX", [x for q in row for x in (a, q)]) for a, row in zip(xa, x_checks))
            s += line("H", xa) + line("MR", za + xa)
        return s

    c = line("R", data + za + xa) + body() + "".join(f"DETECTOR rec[-{k}]\n" for k in (8, 7, 6, 5))
    c += f"REPEAT {rounds - 1} {{\n" + body()
    c += "".join(f"DETECTOR rec[-{k}] rec[-{k + 8}]\n" for k in range(8, 0, -1)) + "}\n"
    c += line("M", data)  # data 0..8: rec[-9]..rec[-1]; the last Z-ancilla results: rec[-17]..rec[-14]
    for a_rec, row in zip((17, 16, 15, 14), z_checks):
        c += "DETECTOR " + " ".join(f"rec[-{9 - q}]" for q in row) + f" rec[-{a_rec}]\n"
    c += "OBSERVABLE_INCLUDE(0) rec[-9] rec[-8] rec[-7]\n"  # logical Z on the data qubits 0, 1, 2
    return stim.Circuit(c)


def test_classify_qubits_checks_measured_in_hadamard_frame():
    # E.2.10: single-qubit gates on data qubits act on the Paulis of all checks, including those already measured,
    # so that all checks are expressed in the same frame
    if not dist_m4ri._HAS_STIM:
        pytest.skip("stim is not installed")
    for x_in_h_frame, meas_before_h in ((False, False), (True, False), (True, True)):
        circ = _surface_memory_z(x_in_h_frame, meas_before_h)
        circ.detector_error_model()  # deterministic detectors and observable
        t_res = dist_m4ri.classify_qubits_thorough(circ)
        assert t_res["is_css"] is True and t_res["is_rotated_css"] is False and t_res["basis"] == "Z"
        assert t_res["data_qubits"] == list(range(9))
        assert t_res["z_ancillas"] == [9, 10, 11, 12] and t_res["x_ancillas"] == [13, 14, 15, 16]
        _, stripped, kept = dist_m4ri.strip_minority_detectors(circ, "Z", thorough_res=t_res)
        assert (stripped, kept) == (8, 16)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])



