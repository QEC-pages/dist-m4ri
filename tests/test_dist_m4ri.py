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


def test_dem_distance_file():
    dem_file = os.path.join(EXAMPLES_DIR, "surf_d3.dem")
    dist, d_list = dist_m4ri.compute_dem_distance(dem=dem_file, threads=4)
    assert dist == 3

    # With codewords (method=3 bracketing mode collects discovered cws)
    dist, d_list, cws = dist_m4ri.compute_dem_distance(dem=dem_file, do_cws=True, threads=4)
    assert dist == 3
    assert len(cws) > 0
    assert all(len(cw) == 3 for cw in cws)

    # Exhaustive CC scan (method=2) finds all 128 codewords
    dist, d_list, cws_cc = dist_m4ri.compute_dem_distance(dem=dem_file, method=2, wmax=3, do_cws=True, threads=4)
    assert dist == 3
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
    dist, d_list = dist_m4ri.compute_dem_distance(dem=dem, threads=4)
    assert dist == 3


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

    d = dist_m4ri.compute_classical_distance(H, solver="codedistance", num_steps=100)
    assert d == 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
