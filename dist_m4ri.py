"""
dist_m4ri.py: Python wrapper for the multithreaded dist_m4ri distance calculator.

Provides high-level APIs for computing:
- Classical code distance: compute_classical_distance(...)
- CSS quantum code distance: compute_css_distance(...)
- Detector Error Model (DEM) distance: compute_dem_distance(...)

Supports distance caching, codeword export, and fallback/option for the codedistance library.
Since dist_m4ri natively handles multithreading (via POSIX threads and dynamic bracketing),
no Python-level threading or subprocess-per-core logic is necessary.
"""

import os
import sys
import time
import random
import shutil
import tempfile
import threading
import subprocess
from pathlib import Path
from typing import List, Tuple, Union, Optional, Dict, Any

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import issparse, csr_array, csr_matrix

try:
    import stim
    _HAS_STIM = True
except ImportError:
    _HAS_STIM = False

try:
    import codedistance
    _HAS_CODEDISTANCE = True
except ImportError:
    _HAS_CODEDISTANCE = False

# Global cache for distance results
_distance_cache: Dict[Any, Any] = {}
_use_distance_cache: bool = True

# Aliases for backward compatibility with vecdec.py
_css_distance_cache = _distance_cache
_use_css_distance_cache = _use_distance_cache


def clear_distance_cache():
    """Clears all cached distance calculations."""
    global _distance_cache
    _distance_cache.clear()


def enable_distance_cache():
    """Enables distance caching."""
    global _use_distance_cache
    _use_distance_cache = True


def disable_distance_cache():
    """Disables distance caching."""
    global _use_distance_cache
    _use_distance_cache = False


# Backward-compatibility aliases
clear_css_distance_cache = clear_distance_cache
enable_css_distance_cache = enable_distance_cache
disable_css_distance_cache = disable_distance_cache


def get_sparse_array_state(A) -> Any:
    """Returns a hashable state representing a numpy/scipy matrix for cache keys."""
    if A is None:
        return None
    if isinstance(A, (str, Path)):
        return str(A)
    if isinstance(A, np.ndarray):
        return (A.shape, A.dtype.str, A.tobytes())
    if issparse(A):
        csr = A.tocsr()
        return (csr.shape, csr.data.tobytes(), csr.indices.tobytes(), csr.indptr.tobytes())
    if hasattr(A, 'tobytes'):
        return A.tobytes()
    return str(A)


def create_unique_file(directory: Union[str, Path] = "tmp", extension: str = ".tmp") -> str:
    """Creates a unique temporary file path and ensures the parent directory exists."""
    os.makedirs(directory, exist_ok=True)
    fd, path = tempfile.mkstemp(suffix=extension, dir=directory)
    os.close(fd)
    return path


def read_sparse_vectors(filepath: str) -> List[List[int]]:
    """
    Reads a list of sparse vectors from a text file in NZLIST format,
    converting from 1-based indexing (in the file) to 0-based indexing (in Python).

    Args:
        filepath (str): The path to the text file.

    Returns:
        list of list of int: A list where each element is a 0-based sparse vector.
    """
    sparse_vectors = []
    if not os.path.exists(filepath):
        return sparse_vectors

    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        if first_line != '%% NZLIST':
            raise ValueError(f"Invalid file format in {filepath}: Missing '%% NZLIST' header.")

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            try:
                parts = list(map(int, line.split()))
            except ValueError:
                raise ValueError(f"Non-integer data found on line {line_num}: {line}")

            stated_length = parts[0]
            vector_elements = [x - 1 for x in parts[1:]]
            if len(vector_elements) != stated_length:
                raise ValueError(
                    f"Length mismatch on line {line_num}. "
                    f"Expected {stated_length} elements, but found {len(vector_elements)}."
                )
            sparse_vectors.append(vector_elements)

    return sparse_vectors


def find_dist_m4ri_binary(custom_path: Optional[str] = None) -> str:
    """Finds the dist_m4ri executable."""
    if custom_path and os.path.isfile(custom_path) and os.access(custom_path, os.X_OK):
        return os.path.abspath(custom_path)

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(pkg_dir, "src", "dist_m4ri"),
        os.path.join(pkg_dir, "dist_m4ri"),
        os.path.join(pkg_dir, "bin", "dist_m4ri"),
        os.path.join(pkg_dir, "..", "dist-m4ri", "src", "dist_m4ri"),
        os.path.join(os.getcwd(), "src", "dist_m4ri"),
        os.path.join(os.getcwd(), "dist_m4ri"),
        os.path.join(os.getcwd(), "bin", "dist_m4ri"),
    ]

    for cand in candidates:
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return os.path.abspath(cand)

    which_path = shutil.which("dist_m4ri")
    if which_path:
        return which_path

    raise FileNotFoundError(
        "Could not find executable 'dist_m4ri'. Please run 'make -C src' to build it."
    )


def parse_dist_m4ri_output(stdout: str) -> Tuple[int, int, int]:
    """
    Parses the standard output of dist_m4ri.
    Expected format on stdout: "dmin dmax rw_steps", "dmin dmax", or a single integer.
    
    Returns:
        tuple (dmin, dmax, rw_steps)
    """
    lines = stdout.strip().split('\n')
    for line in reversed(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                return int(parts[0]), int(parts[1]), int(parts[2])
            except ValueError:
                continue
        elif len(parts) == 2:
            try:
                return int(parts[0]), int(parts[1]), 0
            except ValueError:
                continue
        elif len(parts) == 1:
            try:
                val = int(parts[0])
                return val, val, 0
            except ValueError:
                continue

    raise RuntimeError(f"Could not parse dist_m4ri output: {stdout}")


def run_dist_m4ri(
    dist_m4ri_path: Optional[str] = None,
    method: int = 3,
    finH: Optional[str] = None,
    finG: Optional[str] = None,
    finL: Optional[str] = None,
    fdem: Optional[str] = None,
    dmin: int = 0,
    dmax: int = 0,
    wmax: int = 0,
    wmin: int = 1,
    dexp: int = 0,
    dest: int = 0,
    steps: Optional[int] = None,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    smax: Optional[int] = None,
    noscan: int = 0,
    classical: int = -1,
    dW: int = -1,
    maxC: int = 0,
    pmin: float = 0.0,
    outC: Optional[str] = None,
    seed: int = 0,
    debug: int = 0,
    stop_event: Optional[threading.Event] = None
) -> Tuple[int, int, int]:
    """
    Low-level invocation of the multithreaded dist_m4ri binary.
    
    Returns:
        tuple (dmin, dmax, rw_steps)
    """
    exec_path = find_dist_m4ri_binary(dist_m4ri_path)
    cmd = [exec_path, f"debug={debug}", f"method={method}"]

    if finH: cmd.append(f"finH={finH}")
    if finG: cmd.append(f"finG={finG}")
    if finL: cmd.append(f"finL={finL}")
    if fdem: cmd.append(f"fdem={fdem}")
    if dmin > 0: cmd.append(f"dmin={dmin}")
    if dmax > 0: cmd.append(f"dmax={dmax}")
    if wmax > 0: cmd.append(f"wmax={wmax}")
    if wmin > 1: cmd.append(f"wmin={wmin}")
    if dexp > 0: cmd.append(f"dexp={dexp}")
    elif dest > 0: cmd.append(f"dest={dest}")
    if steps is not None and steps > 0: cmd.append(f"steps={steps}")
    if threads is not None and threads > 0: cmd.append(f"threads={threads}")
    if timeout > 0: cmd.append(f"timeout={timeout}")
    if smax is not None: cmd.append(f"smax={smax}")
    if noscan: cmd.append(f"noscan={noscan}")
    if classical >= 0: cmd.append(f"classical={classical}")
    if dW >= 0: cmd.append(f"dW={dW}")
    if maxC > 0: cmd.append(f"maxC={maxC}")
    if pmin > 0.0: cmd.append(f"pmin={pmin}")
    if outC: cmd.append(f"outC={outC}")
    if seed != 0: cmd.append(f"seed={seed}")

    if debug & 2:
        print(f"[dist_m4ri] Running: {' '.join(cmd)}")

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if stop_event is not None:
        while proc.poll() is None:
            if stop_event.is_set():
                proc.terminate()
                try:
                    proc.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    proc.kill()
                raise RuntimeError("dist_m4ri execution cancelled by stop_event")
            time.sleep(0.05)
        stdout, stderr = proc.communicate()
    else:
        stdout, stderr = proc.communicate()

    if proc.returncode != 0:
        raise RuntimeError(f"dist_m4ri failed with exit code {proc.returncode}:\n{stderr}")

    return parse_dist_m4ri_output(stdout)


def _matrix_to_file(matrix, extension: str = ".mtx", temp_dir: str = "tmp") -> str:
    """Helper to convert a matrix (numpy or scipy sparse) or file path to an MTX file path."""
    if isinstance(matrix, (str, Path)):
        return str(matrix)

    path = create_unique_file(directory=temp_dir, extension=extension)
    if issparse(matrix):
        csr = matrix.astype(np.int8)
        mmwrite(path, csr, symmetry='general')
    else:
        mat_arr = np.asarray(matrix, dtype=np.int8)
        csr = csr_matrix(mat_arr)
        mmwrite(path, csr, symmetry='general')
    return path


def compute_classical_distance(
    H: Any,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmax: int = 0,
    dW: int = -1,
    maxC: int = 0,
    do_cws: bool = False,
    solver: str = "dist_m4ri",
    codedistance_method: str = "QDistEvol",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0
) -> Union[int, Tuple[int, List[List[int]]]]:
    """
    Computes the minimum distance of a classical linear code given parity check matrix H.

    Args:
        H: Parity check matrix (numpy array, scipy sparse matrix, or file path).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance.
        d_max / dmax: Known upper bound on distance.
        wmax: Maximum weight to search in CC.
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        do_cws: Whether to return extracted codewords.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug level flags.

    Returns:
        dist or (dist, cws) if do_cws is True
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max

    global _distance_cache, _use_distance_cache
    cache_key = None
    if _use_distance_cache:
        try:
            h_state = get_sparse_array_state(H)
            cache_key = ("classical", h_state, method, d_exp, eff_dmin, eff_dmax, wmax, dW, maxC, do_cws, solver, codedistance_method)
            if cache_key in _distance_cache:
                if debug & 4:
                    print("[dist_m4ri] Cache hit for classical distance!")
                return _distance_cache[cache_key]
        except Exception:
            cache_key = None

    if solver == "codedistance":
        if not _HAS_CODEDISTANCE:
            raise ImportError("codedistance library requested but not installed.")
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        H_mat = H.toarray() if issparse(H) else (np.asarray(H, dtype=np.int8) if isinstance(H, np.ndarray) else None)
        res = codedistance.codeDistance(
            H_mat, None, tB=1, method=codedistance_method, params=params,
            seed=seed if seed != 0 else None
        )
        d = res.get("d", -1)
        cws = []
        if do_cws and res.get("L") is not None and np.sum(res["L"]) > 0:
            cws = [list(np.where(res["L"])[0])]
        res_val = (d, cws) if do_cws else d
        if _use_distance_cache and cache_key is not None:
            _distance_cache[cache_key] = res_val
        return res_val

    # Use native multithreaded dist_m4ri
    temp_files = []
    try:
        if isinstance(H, (str, Path)) and os.path.exists(str(H)):
            file_H = str(H)
        else:
            file_H = _matrix_to_file(H, extension="_H.mtx")
            temp_files.append(file_H)

        outC = None
        if do_cws:
            outC = create_unique_file(extension="_cws.nz")
            temp_files.append(outC)

        dmin, dmax, _ = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            finH=file_H,
            classical=1,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmax=wmax,
            dexp=d_exp,
            steps=num_steps,
            threads=threads,
            timeout=timeout,
            dW=dW,
            maxC=maxC,
            outC=outC,
            seed=seed,
            debug=debug
        )

        dist = dmin if (dmin == dmax or dmax == 0) else dmax
        cws = []
        if do_cws and outC and os.path.exists(outC):
            cws = read_sparse_vectors(outC)
            cws.sort(key=len)

        res_val = (dist, cws) if do_cws else dist
        if _use_distance_cache and cache_key is not None:
            _distance_cache[cache_key] = res_val
        return res_val

    finally:
        for f in temp_files:
            if os.path.exists(f):
                try: os.remove(f)
                except OSError: pass


def compute_css_distance(
    Hx: Any,
    Hz: Any,
    Lx: Optional[Any] = None,
    Lz: Optional[Any] = None,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmax: int = 0,
    dW: int = -1,
    maxC: int = 0,
    do_cws: bool = False,
    solver: str = "dist_m4ri",
    codedistance_method: str = "QDistEvol",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0,
    **kwargs
) -> Tuple[Any, ...]:
    """
    Computes CSS quantum code distance d = min(d_X, d_Z).

    Args:
        Hx: X-stabilizer parity check matrix.
        Hz: Z-stabilizer parity check matrix.
        Lx: Optional X-logical operator matrix (alternative to Hz as finG).
        Lz: Optional Z-logical operator matrix (alternative to Hx as finG).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance, inclusive.
        d_max / dmax: Known upper bound on distance, inclusive.
        wmax: Maximum weight to search in CC.
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        do_cws: Whether to return extracted X and Z codewords.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug level flags.

    Returns:
        tuple (dist, dist_list_X, dist_list_Z, cws_X, cws_Z) if do_cws
        else (dist, dist_list_X, dist_list_Z)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max

    can_compute_Z = Hx is not None and (hasattr(Hx, 'shape') and Hx.shape[0] > 0 if not isinstance(Hx, str) else True)
    can_compute_X = Hz is not None and (hasattr(Hz, 'shape') and Hz.shape[0] > 0 if not isinstance(Hz, str) else True)

    if not can_compute_Z and not can_compute_X:
        raise ValueError("Cannot compute CSS distance: Both Hx and Hz are empty.")

    global _distance_cache, _use_distance_cache
    cache_key = None
    if _use_distance_cache:
        try:
            hx_state = get_sparse_array_state(Hx) if can_compute_Z else None
            hz_state = get_sparse_array_state(Hz) if can_compute_X else None
            cache_key = (
                "css", hx_state, hz_state, method, d_exp, eff_dmin, eff_dmax, wmax, dW, maxC,
                do_cws, solver, codedistance_method
            )
            if cache_key in _distance_cache:
                if debug & 4:
                    print("[dist_m4ri] Cache hit for CSS distance!")
                return _distance_cache[cache_key]
        except Exception:
            cache_key = None

    if solver == "codedistance":
        if not _HAS_CODEDISTANCE:
            raise ImportError("codedistance library requested but not installed.")
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        cws_Z, cws_X = [], []
        dist_Z, dist_X = None, None
        dist_list_Z, dist_list_X = [], []

        Hx_mat = Hx.toarray() if issparse(Hx) else (np.asarray(Hx, dtype=np.int8) if isinstance(Hx, np.ndarray) else None)
        Hz_mat = Hz.toarray() if issparse(Hz) else (np.asarray(Hz, dtype=np.int8) if isinstance(Hz, np.ndarray) else None)

        if can_compute_Z:
            res_Z = codedistance.CSScodeDistance(
                Hx_mat, Hz_mat, method=codedistance_method, params=params.copy(),
                component="Z", seed=seed if seed != 0 else None
            )
            dist_Z = res_Z.get("d", -1)
            dist_list_Z = [dist_Z] if dist_Z > 0 else []
            if do_cws and res_Z.get("L") is not None and np.sum(res_Z["L"]) > 0:
                cws_Z = [list(np.where(res_Z["L"])[0])]

        if can_compute_X:
            res_X = codedistance.CSScodeDistance(
                Hx_mat, Hz_mat, method=codedistance_method, params=params.copy(),
                component="X", seed=seed if seed != 0 else None
            )
            dist_X = res_X.get("d", -1)
            dist_list_X = [dist_X] if dist_X > 0 else []
            if do_cws and res_X.get("L") is not None and np.sum(res_X["L"]) > 0:
                cws_X = [list(np.where(res_X["L"])[0])]

        if dist_X is not None and dist_Z is not None:
            dist = min(dist_Z, dist_X) if (dist_Z > 0 and dist_X > 0) else max(dist_Z, dist_X)
        elif dist_X is not None:
            dist = dist_X
        else:
            dist = dist_Z

        res_tuple = (dist, dist_list_X, dist_list_Z, cws_X, cws_Z) if do_cws else (dist, dist_list_X, dist_list_Z)
        if _use_distance_cache and cache_key is not None:
            _distance_cache[cache_key] = res_tuple
        return res_tuple

    # Native multithreaded dist_m4ri
    temp_files = []
    try:
        file_Hx = _matrix_to_file(Hx, extension="_Hx.mtx") if can_compute_Z else None
        file_Hz = _matrix_to_file(Hz, extension="_Hz.mtx") if can_compute_X else None
        file_Lx = _matrix_to_file(Lx, extension="_Lx.mtx") if Lx is not None else None
        file_Lz = _matrix_to_file(Lz, extension="_Lz.mtx") if Lz is not None else None

        for f in (file_Hx, file_Hz, file_Lx, file_Lz):
            if f and not isinstance(f, (str, Path)) or (f and not os.path.exists(f)):
                pass
            elif f and f.startswith(tempfile.gettempdir()):
                temp_files.append(f)

        outZ = create_unique_file(extension="_Z.nz") if (do_cws and can_compute_Z) else None
        outX = create_unique_file(extension="_X.nz") if (do_cws and can_compute_X) else None
        if outZ: temp_files.append(outZ)
        if outX: temp_files.append(outX)

        dist_Z, dist_X = None, None
        dist_list_Z, dist_list_X = [], []
        cws_Z, cws_X = [], []

        # Z-distance: Hx as finH, Hz as finG (or Lz as finL)
        if can_compute_Z:
            dmin_z, dmax_z, _ = run_dist_m4ri(
                dist_m4ri_path=dist_m4ri,
                method=method,
                finH=file_Hx,
                finG=file_Hz if file_Lz is None else None,
                finL=file_Lz,
                dmin=eff_dmin,
                dmax=eff_dmax,
                wmax=wmax,
                dexp=d_exp,
                steps=num_steps,
                threads=threads,
                timeout=timeout,
                dW=dW,
                maxC=maxC,
                outC=outZ,
                seed=seed,
                debug=debug
            )
            dist_Z = dmin_z if (dmin_z == dmax_z or dmax_z == 0) else dmax_z
            dist_list_Z = [dist_Z] if dist_Z > 0 else []
            if do_cws and outZ and os.path.exists(outZ):
                cws_Z = read_sparse_vectors(outZ)
                cws_Z.sort(key=len)

        # X-distance: Hz as finH, Hx as finG (or Lx as finL)
        if can_compute_X:
            dmin_x, dmax_x, _ = run_dist_m4ri(
                dist_m4ri_path=dist_m4ri,
                method=method,
                finH=file_Hz,
                finG=file_Hx if file_Lx is None else None,
                finL=file_Lx,
                dmin=eff_dmin,
                dmax=eff_dmax,
                wmax=wmax,
                dexp=d_exp,
                steps=num_steps,
                threads=threads,
                timeout=timeout,
                dW=dW,
                maxC=maxC,
                outC=outX,
                seed=seed,
                debug=debug
            )
            dist_X = dmin_x if (dmin_x == dmax_x or dmax_x == 0) else dmax_x
            dist_list_X = [dist_X] if dist_X > 0 else []
            if do_cws and outX and os.path.exists(outX):
                cws_X = read_sparse_vectors(outX)
                cws_X.sort(key=len)

        if dist_X is not None and dist_Z is not None:
            dist = min(dist_Z, dist_X) if (dist_Z > 0 and dist_X > 0) else max(dist_Z, dist_X)
        elif dist_X is not None:
            dist = dist_X
        else:
            dist = dist_Z

        res_tuple = (dist, dist_list_X, dist_list_Z, cws_X, cws_Z) if do_cws else (dist, dist_list_X, dist_list_Z)
        if _use_distance_cache and cache_key is not None:
            _distance_cache[cache_key] = res_tuple
        return res_tuple

    finally:
        for f in temp_files:
            if f and os.path.exists(f):
                try: os.remove(f)
                except OSError: pass


def compute_dem_distance(
    dem: Optional[Any] = None,
    circuit: Optional[Any] = None,
    dist_m4ri: Optional[str] = None,
    method: int = 3,
    threads: Optional[int] = None,
    timeout: float = 60.0,
    num_steps: Optional[int] = None,
    d_exp: int = 0,
    d_min: int = 0,
    d_max: int = 0,
    dmin: int = 0,
    dmax: int = 0,
    wmax: int = 0,
    dW: int = -1,
    maxC: int = 0,
    pmin: float = 0.0,
    noscan: int = 0,
    do_cws: bool = False,
    solver: str = "dist_m4ri",
    codedistance_method: str = "UndetectableErrorStim",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0,
    **kwargs
) -> Tuple[Any, ...]:
    """
    Computes minimum graph/hypergraph distance of a Stim Detector Error Model (DEM).

    Args:
        dem: stim.DetectorErrorModel object or path to .dem file.
        circuit: stim.Circuit object (converted to DEM).
        dist_m4ri: Path to dist_m4ri executable (optional).
        method: Solver method (1=RW, 2=CC, 3=Bracketing default).
        threads: Number of worker threads.
        timeout: Execution timeout in seconds.
        num_steps: Maximum RW steps.
        d_exp: Expected distance estimate.
        d_min / dmin: Known lower bound on distance, inclusive.
        d_max / dmax: Known upper bound on distance, inclusive.
        wmax: Maximum weight to search in CC.
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        pmin: Probability cutoff for error mechanisms in DEM.
        noscan: Skip CC scan loop if 1.
        do_cws: Whether to return extracted error mechanisms / codewords.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug level flags.

    Returns:
        tuple (dist, dist_list, cws) if do_cws else (dist, dist_list)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max

    if dem is None and circuit is not None:
        if hasattr(circuit, 'detector_error_model'):
            dem = circuit.detector_error_model(decompose_errors=True)
        else:
            raise ValueError("Provided circuit object does not have detector_error_model() method.")

    if dem is None:
        raise ValueError("Either 'dem' or 'circuit' must be provided.")

    if solver == "codedistance":
        if not _HAS_CODEDISTANCE:
            raise ImportError("codedistance library requested but not installed.")
        params = dict(codedistance_params or {})
        params.setdefault("filterCircuit", False)
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        if circuit is not None:
            res = codedistance.circuitDistance(
                circuit, method=codedistance_method, params=params,
                seed=seed if seed != 0 else None
            )
        else:
            H, L, priors = codedistance.StimDEM2HL(dem)
            if "priors" not in params and len(priors) > 0:
                params["priors"] = priors
            res = codedistance.codeDistance(
                H, L, tB=1, method=codedistance_method, params=params,
                seed=seed if seed != 0 else None
            )
        d = res.get("d", -1)
        dist_list = [d] if d > 0 else []
        if do_cws:
            cws = [list(np.where(res["L"])[0])] if (res.get("L") is not None and np.sum(res["L"]) > 0) else []
            return d, dist_list, cws
        return d, dist_list

    # Native multithreaded dist_m4ri directly on DEM
    temp_files = []
    try:
        if isinstance(dem, (str, Path)) and os.path.exists(str(dem)):
            file_dem = str(dem)
        else:
            file_dem = create_unique_file(extension=".dem")
            temp_files.append(file_dem)
            if hasattr(dem, 'flattened'):
                dem.flattened().to_file(file_dem)
            elif hasattr(dem, 'to_file'):
                dem.to_file(file_dem)
            else:
                with open(file_dem, 'w') as f:
                    f.write(str(dem))

        outC = None
        if do_cws:
            outC = create_unique_file(extension="_out.nz")
            temp_files.append(outC)

        dmin, dmax, _ = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            fdem=file_dem,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmax=wmax,
            dexp=d_exp,
            steps=num_steps,
            threads=threads,
            timeout=timeout,
            smax=0,
            noscan=noscan,
            pmin=pmin,
            dW=dW,
            maxC=maxC,
            outC=outC,
            seed=seed,
            debug=debug
        )

        dist = dmin if (dmin == dmax or dmax == 0) else dmax
        dist_list = [dist] if dist > 0 else []
        cws = []
        if do_cws and outC and os.path.exists(outC):
            cws = read_sparse_vectors(outC)
            cws.sort(key=len)

        if do_cws:
            return dist, dist_list, cws
        return dist, dist_list

    finally:
        for f in temp_files:
            if f and os.path.exists(f):
                try: os.remove(f)
                except OSError: pass
