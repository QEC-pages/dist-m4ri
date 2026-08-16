"""
Unit tests for dist_m4ri.py Python wrapper.
"""

import os
import sys
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

    # Exhaustive CC scan (method=2) finds all 128 codewords
    dist, d_info, cws_cc = dist_m4ri.compute_dem_distance(dem=dem_file, method=2, wmax=3, do_cws=True, threads=4)
    assert dist == 3
    assert d_info == [3, 3, 0]
    assert len(cws_cc) == 128
    assert all(len(cw) == 3 for cw in cws_cc)


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

    # Method 1 (RW): rw_steps reported
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    dmin, dmax, rw_steps = dist_m4ri.run_dist_m4ri(method=1, fdem=dem_file, steps=100, threads=4)
    assert dmin == 1
    assert dmax == 3
    assert rw_steps >= 100


def test_dmin_dmax_parameters():
    # Test dmin/dmax in run_dist_m4ri
    h_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    l_file = os.path.join(EXAMPLES_DIR, "surf_d5_L.mmx")
    dmin, dmax, rw_steps = dist_m4ri.run_dist_m4ri(method=3, finH=h_file, finL=l_file, dmin=4, dmax=5, timeout=5, threads=4)
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
    assert len(data1) == 1
    key = list(data1.keys())[0]
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
    d_css, dx_info, dz_info = dist_m4ri.compute_css_distance(
        Hx=hx_file, Hz=hx_file, Lz=hz_file, Lx=hz_file,
        d_exp=5, threads=4, cache_file=json_file
    )
    assert d_css == 5
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


def test_wmin_greater_than_dmax():
    hx_file = os.path.join(EXAMPLES_DIR, "surf_d5_H.mmx")
    # When dmax=5 and wmin=7, wmin >= dmax should terminate immediately without search
    d, d_info = dist_m4ri.compute_classical_distance(
        H=hx_file, dmax=5, wmin=7, return_info=True, threads=4
    )
    assert d == 5
    assert d_info == [0, 5, 0]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
