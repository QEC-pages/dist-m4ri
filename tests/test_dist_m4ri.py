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
    assert args1["start"] == 2

    # Auto-infer classical = 1 when only finH is given
    args2 = dist_m4ri.parse_cli_args(["finH=h.mtx", "method=2"])
    assert args2["classical"] == 1
    assert args2["finH"] == "h.mtx"
    assert args2["finG"] is None

    # Auto-infer classical = 0 when fdem is given
    args3 = dist_m4ri.parse_cli_args(["fdem=model.dem"])
    assert args3["classical"] == 0


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
    assert "0.9.0" in captured.out


def test_cli_binary_compatibility_silent(capsys):
    ret = dist_m4ri.main(["--help"])
    assert ret == 0
    captured = capsys.readouterr()
    # When binary is found and up to date, stderr should be silent (no warnings)
    assert "Warning:" not in captured.err
    assert "0.9.0" in captured.out


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
    assert "expected >= 0.9.0" in older_warn


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
        d, d_info = dist_m4ri.compute_dem_distance(
            dem=dem_file, method=1, num_steps=2000, ksub=32, kwin=48,
            win_mode=0, min_hits=3, cov_cws=2, refresh=50, threads=4, seed=42
        )
        assert d == 3
        assert d_info[1] == 3
        assert 0 < d_info[2] < 2000

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
        # Default for CSS .stim is --simple
        ret = dist_m4ri.main([str(stim_file), "-v", "--no-cache", "threads=2"])
        assert ret == 0
        out_simple = capsys.readouterr().out
        assert "mode: simple" in out_simple
        assert "kept=16, stripped=8" in out_simple

        # Explicit --full keeps all 24 detectors
        ret = dist_m4ri.main(
            [str(stim_file), "--full", "-v", "--no-cache", "threads=2"]
        )
        assert ret == 0
        out_full = capsys.readouterr().out
        assert "mode: full" in out_full
        assert "kept=24, stripped=0" in out_full
    finally:
        dist_m4ri.enable_distance_cache()


def test_stim_xzzx_rotated_css(tmp_path):
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
        d_simp, info_simp = dist_m4ri.compute_dem_distance(
            dem=str(stim_path), simple=True, threads=2
        )
        assert d_simp == 3
        assert info_simp[0] == 3 and info_simp[1] == 3

        dist_m4ri.clear_distance_cache()
        d_full, info_full = dist_m4ri.compute_dem_distance(
            dem=str(stim_path), full=True, threads=2
        )
        assert d_full == 3
        assert info_full[0] == 3 and info_full[1] == 3


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

    # CLI argument parsing checks for --rounds, --simple, --full
    args1 = dist_m4ri.parse_cli_args(["--rounds", str(stim_with_repeat)])
    assert args1["rounds"] == 2
    assert args1["fdem"] == str(stim_with_repeat)

    args2 = dist_m4ri.parse_cli_args(["--rounds", "4", "--simple", str(stim_with_repeat)])
    assert args2["rounds"] == 4
    assert args2["simple"] is True
    assert args2["full"] is False
    assert args2["fdem"] == str(stim_with_repeat)

    args3 = dist_m4ri.parse_cli_args(["rounds=3", "--full", "basis=X", str(stim_with_repeat)])
    assert args3["rounds"] == 3
    assert args3["full"] is True
    assert args3["simple"] is False
    assert args3["basis"] == "X"

    dist_m4ri.clear_distance_cache()
    dist_m4ri.disable_distance_cache()
    try:
        # Bare --rounds sets rounds=2 (REPEAT 2 block -> 16 Z-detectors kept, 8 X-detectors stripped)
        ret = dist_m4ri.main(
            [str(stim_with_repeat), "--rounds", "-v", "--no-cache", "threads=2"]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "rounds=2" in out
        assert "kept=16, stripped=8" in out

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
        # 2. Bare --out-dem and --out-stim (default --simple, default out-dir = input file dir)
        ret = dist_m4ri.main(
            [str(stim_file), "--rounds", "2", "--out-dem", "--out-stim", "--no-cache", "threads=2"]
        )
        assert ret == 0
        exp_simp_dem = in_dir / "surf_d3_Z_simp.dem"
        exp_simp_stim = in_dir / "surf_d3_Z_simp.stim"
        assert exp_simp_dem.exists()
        assert exp_simp_stim.exists()

        saved_simp_circ = stim.Circuit.from_file(str(exp_simp_stim))
        assert dist_m4ri.has_noise(saved_simp_circ)
        assert saved_simp_circ.num_detectors == 16
        saved_simp_dem = stim.DetectorErrorModel.from_file(str(exp_simp_dem))
        assert saved_simp_dem.num_detectors == 16

        # 3. Bare --out-dem and --out-stim with --full and explicit --out-dir
        custom_dir = tmp_path / "custom_out"
        ret_full = dist_m4ri.main(
            [
                str(stim_file),
                "--full",
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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])



