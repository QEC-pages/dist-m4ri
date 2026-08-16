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
import json
import time
import random
import shutil
import hashlib
import tempfile
import threading
import subprocess
from pathlib import Path
from typing import List, Tuple, Union, Optional, Dict, Any

_codedistance_mod = None
_stim_mod = None


def _get_codedistance():
    """Lazily imports the codedistance library only when requested."""
    global _codedistance_mod
    if _codedistance_mod is None:
        try:
            import codedistance
            _codedistance_mod = codedistance
        except ImportError:
            raise ImportError("codedistance library is requested but not installed.")
    return _codedistance_mod


def _get_stim():
    """Lazily imports stim only when requested."""
    global _stim_mod
    if _stim_mod is None:
        try:
            import stim
            _stim_mod = stim
        except ImportError:
            raise ImportError("stim library is requested but not installed.")
    return _stim_mod


def __getattr__(name: str) -> Any:
    """Lazily resolves module attributes without loading heavy dependencies at startup."""
    if name == "_HAS_STIM":
        try:
            import stim
            return True
        except ImportError:
            return False
    if name == "_HAS_CODEDISTANCE":
        try:
            import codedistance
            return True
        except ImportError:
            return False
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

# Global cache for distance results
_distance_cache: Dict[str, Any] = {}
_use_distance_cache: bool = True
_distance_cache_file: Optional[str] = None

# Aliases for backward compatibility with vecdec.py
_css_distance_cache = _distance_cache
_use_css_distance_cache = _use_distance_cache


def set_distance_cache_file(filepath: Optional[Union[str, Path]] = None) -> None:
    """
    Sets the default JSON file for persistent distance caching.
    If the file exists, its contents are loaded into memory.
    """
    global _distance_cache_file
    if filepath is not None:
        _distance_cache_file = str(Path(filepath).resolve())
        load_distance_cache(_distance_cache_file)
    else:
        _distance_cache_file = None


def load_distance_cache(filepath: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """
    Loads distance cache from a JSON file into memory.
    """
    global _distance_cache, _distance_cache_file
    target_file = str(Path(filepath).resolve()) if filepath is not None else _distance_cache_file
    if target_file and os.path.isfile(target_file):
        try:
            with open(target_file, "r") as f:
                data = json.load(f)
            if isinstance(data, dict):
                _distance_cache.update(data)
        except Exception as e:
            sys.stderr.write(f"# Warning: Failed to load distance cache from {target_file}: {e}\n")
    return _distance_cache


def save_distance_cache(filepath: Optional[Union[str, Path]] = None) -> None:
    """
    Saves the in-memory distance cache to a JSON file.
    Uses atomic write via a temporary file to prevent corruption.
    """
    global _distance_cache, _distance_cache_file
    target_file = str(Path(filepath).resolve()) if filepath is not None else _distance_cache_file
    if not target_file:
        return

    parent_dir = os.path.dirname(os.path.abspath(target_file)) or "."
    os.makedirs(parent_dir, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(suffix=".tmp", prefix="dist_cache_", dir=parent_dir)
    try:
        with open(fd, "w") as f:
            json.dump(_distance_cache, f, indent=2)
        os.replace(temp_path, target_file)
    except Exception as e:
        if os.path.exists(temp_path):
            try: os.remove(temp_path)
            except OSError: pass
        sys.stderr.write(f"# Warning: Failed to save distance cache to {target_file}: {e}\n")


def clear_distance_cache(cache_file: Optional[Union[str, Path]] = None, clear_file: bool = False) -> None:
    """
    Clears all cached distance calculations from memory, and optionally deletes the persistent JSON file.
    """
    global _distance_cache, _distance_cache_file
    _distance_cache.clear()
    target_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    if clear_file and target_file and os.path.isfile(target_file):
        try:
            os.remove(target_file)
        except OSError:
            pass


def enable_distance_cache() -> None:
    """Enables distance caching."""
    global _use_distance_cache
    _use_distance_cache = True


def disable_distance_cache() -> None:
    """Disables distance caching."""
    global _use_distance_cache
    _use_distance_cache = False


def get_distance_cache() -> Dict[str, Any]:
    """Returns the global distance cache dictionary."""
    global _distance_cache
    return _distance_cache


def format_bounds_list(dmin: int, dmax: int, num_rw: int) -> List[int]:
    """
    Returns [dmin, dmax, num_rw] according to:
    - [d, d, 0] if known exactly
    - [dmin, 0, 0] if there is no upper bound (dmax == 0)
    - [0, dmax, num_rw] if there is no lower bound (dmin <= 1)
    - [dmin, dmax, num_rw] otherwise
    """
    eff_dmin = dmin if dmin > 1 else 0
    eff_dmax = dmax if dmax > 0 else 0
    eff_rw = num_rw if num_rw > 0 else 0
    if eff_dmin > 0 and eff_dmin == eff_dmax:
        return [eff_dmin, eff_dmax, 0]
    elif eff_dmax == 0:
        return [eff_dmin, 0, 0]
    elif eff_dmin == 0:
        return [0, eff_dmax, eff_rw]
    else:
        return [eff_dmin, eff_dmax, eff_rw]


def format_bounds_str(bounds: List[int]) -> str:
    """Formats a bounds list [dmin, dmax, num_rw] as a string, stating '(exact)' if bounds coincide."""
    dmin, dmax, num_rw = bounds[0], bounds[1], bounds[2]
    if dmin > 0 and dmin == dmax:
        return f"{dmin} {dmax} {num_rw} (exact)"
    return f"{dmin} {dmax} {num_rw}"


def get_cached_distance(
    H: Optional[Any] = None,
    Hx: Optional[Any] = None,
    Hz: Optional[Any] = None,
    Lx: Optional[Any] = None,
    Lz: Optional[Any] = None,
    dem: Optional[Any] = None,
    circuit: Optional[Any] = None,
    pmin: float = 0.0,
    cache_file: Optional[Union[str, Path]] = None
) -> Optional[Dict[str, Any]]:
    """
    Retrieves the cached distance entry (including bounds and cumulative rw_steps)
    for a given code matrix, CSS code, or DEM.

    Returns:
        dict with keys {"dist", "dmin", "dmax", "rw_steps", ...} or None if not cached.
    """
    global _distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    if eff_cache_file:
        load_distance_cache(eff_cache_file)

    if H is not None:
        key = f"classical:{get_sparse_array_state(H)}"
        entry = _distance_cache.get(key)
        if entry:
            entry = dict(entry)
            entry["d_info"] = format_bounds_list(entry.get("dmin", 0), entry.get("dmax", 0), entry.get("rw_steps", 0))
        return entry
    elif Hx is not None or Hz is not None:
        hx_st = get_sparse_array_state(Hx) if Hx is not None else "none"
        hz_st = get_sparse_array_state(Hz) if Hz is not None else "none"
        key = f"css:X={hx_st}:Z={hz_st}"
        if Lx is not None or Lz is not None:
            lx_st = get_sparse_array_state(Lx) if Lx is not None else "none"
            lz_st = get_sparse_array_state(Lz) if Lz is not None else "none"
            key = f"{key}:Lx={lx_st}:Lz={lz_st}"
        entry = _distance_cache.get(key)
        if entry:
            entry = dict(entry)
            if "dmin_X" in entry:
                entry["dX"] = format_bounds_list(entry.get("dmin_X", 0), entry.get("dmax_X", 0), entry.get("rw_steps_X", 0))
            if "dmin_Z" in entry:
                entry["dZ"] = format_bounds_list(entry.get("dmin_Z", 0), entry.get("dmax_Z", 0), entry.get("rw_steps_Z", 0))
        return entry
    elif dem is not None or circuit is not None:
        if dem is None and circuit is not None:
            if hasattr(circuit, 'detector_error_model'):
                obj = circuit.detector_error_model(decompose_errors=True)
            else:
                obj = circuit
        else:
            obj = dem
        dem_st = get_sparse_array_state(obj)
        key = f"dem:{dem_st}" if pmin <= 0.0 else f"dem:{dem_st}:pmin={pmin}"
        entry = _distance_cache.get(key)
        if entry:
            entry = dict(entry)
            entry["d_info"] = format_bounds_list(entry.get("dmin", 0), entry.get("dmax", 0), entry.get("rw_steps", 0))
        return entry
    return None


# Backward-compatibility aliases
clear_css_distance_cache = clear_distance_cache
enable_css_distance_cache = enable_distance_cache
disable_css_distance_cache = disable_distance_cache


def get_sparse_array_state(A) -> str:
    """Returns a deterministic string representation for JSON-compatible cache keys."""
    if A is None:
        return "none"
    if isinstance(A, (str, Path)):
        path_str = str(Path(A).resolve())
        if os.path.isfile(path_str):
            try:
                with open(path_str, "rb") as f:
                    content_h = hashlib.sha256(f.read()).hexdigest()
                return f"file:{path_str}:{content_h}"
            except Exception:
                return f"file:{path_str}"
        return path_str
    if hasattr(A, 'shape') and hasattr(A, 'dtype') and hasattr(A, 'tobytes'):
        h = hashlib.sha256(A.tobytes()).hexdigest()
        dtype_str = getattr(A.dtype, 'str', str(A.dtype))
        return f"ndarray:{A.shape}:{dtype_str}:{h}"
    if hasattr(A, 'tocsr'):
        csr = A.tocsr()
        h = hashlib.sha256(csr.data.tobytes() + csr.indices.tobytes() + csr.indptr.tobytes()).hexdigest()
        return f"csr:{csr.shape}:{h}"
    if hasattr(A, 'tobytes'):
        h = hashlib.sha256(A.tobytes()).hexdigest()
        return f"bytes:{h}"
    h = hashlib.sha256(str(A).encode('utf-8')).hexdigest()
    return f"str_sha256:{h}"


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
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return sparse_vectors

    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            return sparse_vectors
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


class DistanceResult:
    """
    Structured result for code distance calculations containing:
    - dmin: lower bound on distance
    - dmax: upper bound on distance (minimum non-trivial codeword weight found)
    - rw_steps: cumulative number of completed random window information sets
    - cws: discovered codewords (if requested)
    """
    def __init__(
        self,
        dmin: int,
        dmax: int,
        rw_steps: int = 0,
        cws: Optional[List[List[int]]] = None,
        cws_X: Optional[List[List[int]]] = None,
        cws_Z: Optional[List[List[int]]] = None,
        dmin_X: Optional[int] = None,
        dmax_X: Optional[int] = None,
        rw_steps_X: Optional[int] = None,
        dmin_Z: Optional[int] = None,
        dmax_Z: Optional[int] = None,
        rw_steps_Z: Optional[int] = None,
    ):
        self.dmin = dmin
        self.dmax = dmax
        self.rw_steps = rw_steps
        self.cws = cws
        self.cws_X = cws_X
        self.cws_Z = cws_Z
        self.dmin_X = dmin_X
        self.dmax_X = dmax_X
        self.rw_steps_X = rw_steps_X
        self.dmin_Z = dmin_Z
        self.dmax_Z = dmax_Z
        self.rw_steps_Z = rw_steps_Z

    @property
    def is_exact(self) -> bool:
        return self.dmin > 0 and self.dmin == self.dmax

    @property
    def dist(self) -> int:
        return self.dmin if self.is_exact else (self.dmax if self.dmax > 0 else self.dmin)

    def __int__(self) -> int:
        return self.dist

    def __index__(self) -> int:
        return self.dist

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, DistanceResult):
            return (self.dmin, self.dmax, self.rw_steps) == (other.dmin, other.dmax, other.rw_steps)
        if isinstance(other, (tuple, list)):
            return tuple(self) == tuple(other)
        if isinstance(other, (int, np.integer)):
            return self.dist == other
        return False

    def __iter__(self):
        if self.cws_X is not None or self.cws_Z is not None:
            return iter((self.dmin, self.dmax, self.rw_steps, self.cws_X or [], self.cws_Z or []))
        if self.cws is not None:
            return iter((self.dmin, self.dmax, self.rw_steps, self.cws))
        return iter((self.dmin, self.dmax, self.rw_steps))

    def __getitem__(self, index: int):
        return tuple(self)[index]

    def __len__(self) -> int:
        return len(tuple(self))

    def __str__(self) -> str:
        if self.is_exact:
            return f"{self.dmin} {self.dmax} {self.rw_steps} (exact)"
        return f"{self.dmin} {self.dmax} {self.rw_steps}"

    def __repr__(self) -> str:
        if self.is_exact:
            return f"DistanceResult({self.dmin} {self.dmax} {self.rw_steps} (exact))"
        return f"DistanceResult(dmin={self.dmin}, dmax={self.dmax}, rw_steps={self.rw_steps})"


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

    if method == 2 and wmax <= 0:
        if dmax > 0:
            wmax = dmax
        elif timeout <= 0.0:
            raise ValueError("either parameter wmax>0 or timeout>0 should be specified for CC method=2.")

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

    import numpy as np
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix, issparse

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
    wmin: int = 1,
    wmax: int = 0,
    dW: int = -1,
    maxC: int = 0,
    do_cws: bool = False,
    return_info: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
    solver: str = "dist_m4ri",
    codedistance_method: str = "QDistEvol",
    codedistance_params: Optional[Dict[str, Any]] = None,
    seed: int = 0,
    debug: int = 0
) -> Any:
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
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        do_cws: Whether to return extracted codewords.
        return_info: If True, return (dist, d_info) or (dist, d_info, cws) where d_info is [dmin, dmax, num_rw].
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug level flags.

    Returns:
        dist or (dist, cws) if do_cws is True (or (dist, d_info) / (dist, d_info, cws) if return_info=True)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max

    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None

    if solver == "codedistance":
        if do_cws:
            raise ValueError("Codeword extraction (do_cws=True) is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
        import numpy as np
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        H_mat = H.toarray() if hasattr(H, 'toarray') else (np.asarray(H, dtype=np.int8) if isinstance(H, (np.ndarray, list)) else None)
        res = codedistance.codeDistance(
            H_mat, None, tB=1, method=codedistance_method, params=params,
            seed=seed if seed != 0 else None
        )
        return res.get("d", -1)

    # Solver is native multithreaded dist_m4ri (supports bounds caching and cumulative RW steps)
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            h_state = get_sparse_array_state(H)
            code_key = f"classical:{h_state}"
            cached_entry = _distance_cache.get(code_key)
            if cached_entry is not None:
                # If exact distance is already proven and not asking for more codewords
                if cached_entry.get("dmin", 0) > 0 and cached_entry.get("dmin") == cached_entry.get("dmax"):
                    if not do_cws or (cached_entry.get("cws") and len(cached_entry["cws"]) > 0):
                        if debug & 4:
                            print("[dist_m4ri] Cache hit for classical distance (exact distance known)!")
                        d_info = format_bounds_list(cached_entry.get("dmin", 0), cached_entry.get("dmax", 0), cached_entry.get("rw_steps", 0))
                        if return_info:
                            return (cached_entry["dist"], d_info, cached_entry.get("cws", [])) if do_cws else (cached_entry["dist"], d_info)
                        return (cached_entry["dist"], cached_entry.get("cws", [])) if do_cws else cached_entry["dist"]
                # Use existing cached bounds to accelerate subsequent runs
                if eff_dmax == 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = cached_entry["dmax"]
                elif eff_dmax > 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = min(eff_dmax, cached_entry["dmax"])
                if eff_dmin <= 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = cached_entry["dmin"]
                elif eff_dmin > 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = max(eff_dmin, cached_entry["dmin"])
        except Exception:
            code_key = None
            cached_entry = None

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

        dmin_res, dmax_res, rw_steps = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            finH=file_H,
            classical=1,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmin=wmin,
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

        dist = dmin_res if (dmin_res == dmax_res or dmax_res == 0) else dmax_res
        cws = []
        if do_cws and outC and os.path.exists(outC):
            cws = read_sparse_vectors(outC)
            cws.sort(key=len)

        d_info = format_bounds_list(dmin_res, dmax_res, rw_steps)

        if _use_distance_cache and code_key is not None:
            prev_steps = cached_entry.get("rw_steps", 0) if cached_entry else 0
            prev_dmax = cached_entry.get("dmax", 0) if cached_entry else 0
            prev_dmin = cached_entry.get("dmin", 0) if cached_entry else 0
            prev_cws = list(cached_entry.get("cws", [])) if cached_entry else []

            total_rw_steps = prev_steps + rw_steps
            best_dmax = min(prev_dmax, dmax_res) if (prev_dmax > 0 and dmax_res > 0) else (dmax_res if dmax_res > 0 else prev_dmax)
            best_dmin = max(prev_dmin, dmin_res)

            combined_cws = prev_cws
            if cws:
                existing_set = {tuple(cw) for cw in combined_cws}
                for cw in cws:
                    if tuple(cw) not in existing_set:
                        combined_cws.append(cw)
                        existing_set.add(tuple(cw))
                combined_cws.sort(key=len)

            d_info = format_bounds_list(best_dmin, best_dmax, total_rw_steps)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "d_info": d_info,
                "cws": combined_cws
            }
            if eff_cache_file:
                save_distance_cache(eff_cache_file)

            if return_info:
                return (dist, d_info, combined_cws) if do_cws else (dist, d_info)
            return (dist, combined_cws) if do_cws else dist

        if return_info:
            return (dist, d_info, cws) if do_cws else (dist, d_info)
        return (dist, cws) if do_cws else dist

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
    wmin: int = 1,
    wmax: int = 0,
    dW: int = -1,
    maxC: int = 0,
    do_cws: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
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
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        do_cws: Whether to return extracted X and Z codewords.
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug level flags.

    Returns:
        tuple (dist, dX_info, dZ_info, cws_X, cws_Z) if do_cws
        else (dist, dX_info, dZ_info)
    """
    eff_dmin = dmin if dmin > 0 else d_min
    eff_dmax = dmax if dmax > 0 else d_max

    can_compute_Z = Hx is not None and (hasattr(Hx, 'shape') and Hx.shape[0] > 0 if not isinstance(Hx, str) else True)
    can_compute_X = Hz is not None and (hasattr(Hz, 'shape') and Hz.shape[0] > 0 if not isinstance(Hz, str) else True)

    if not can_compute_Z and not can_compute_X:
        raise ValueError("Cannot compute CSS distance: Both Hx and Hz are empty.")

    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None

    if solver == "codedistance":
        if do_cws:
            raise ValueError("Codeword extraction (do_cws=True) is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
        import numpy as np
        params = dict(codedistance_params or {})
        if num_steps is not None and "iterCount" not in params:
            params["iterCount"] = num_steps

        dist_Z, dist_X = None, None
        dX_info, dZ_info = None, None

        Hx_mat = Hx.toarray() if hasattr(Hx, 'toarray') else (np.asarray(Hx, dtype=np.int8) if isinstance(Hx, (np.ndarray, list)) else None)
        Hz_mat = Hz.toarray() if hasattr(Hz, 'toarray') else (np.asarray(Hz, dtype=np.int8) if isinstance(Hz, (np.ndarray, list)) else None)

        if can_compute_Z:
            res_Z = codedistance.CSScodeDistance(
                Hx_mat, Hz_mat, method=codedistance_method, params=params.copy(),
                component="Z", seed=seed if seed != 0 else None
            )
            dist_Z = res_Z.get("d", -1)
            dZ_info = format_bounds_list(dist_Z, dist_Z, 0) if dist_Z > 0 else [0, 0, 0]

        if can_compute_X:
            res_X = codedistance.CSScodeDistance(
                Hx_mat, Hz_mat, method=codedistance_method, params=params.copy(),
                component="X", seed=seed if seed != 0 else None
            )
            dist_X = res_X.get("d", -1)
            dX_info = format_bounds_list(dist_X, dist_X, 0) if dist_X > 0 else [0, 0, 0]

        if dist_X is not None and dist_Z is not None:
            dist = min(dist_Z, dist_X) if (dist_Z > 0 and dist_X > 0) else max(dist_Z, dist_X)
        elif dist_X is not None:
            dist = dist_X
        else:
            dist = dist_Z

        return (dist, dX_info, dZ_info)

    # Solver is native multithreaded dist_m4ri
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            hx_state = get_sparse_array_state(Hx) if can_compute_Z else "none"
            hz_state = get_sparse_array_state(Hz) if can_compute_X else "none"
            code_key = f"css:X={hx_state}:Z={hz_state}"
            if Lx is not None or Lz is not None:
                lx_state = get_sparse_array_state(Lx) if Lx is not None else "none"
                lz_state = get_sparse_array_state(Lz) if Lz is not None else "none"
                code_key = f"{code_key}:Lx={lx_state}:Lz={lz_state}"
            cached_entry = _distance_cache.get(code_key)
            if cached_entry is not None:
                # If exact distance is already proven and not asking for more codewords
                if cached_entry.get("dmin", 0) > 0 and cached_entry.get("dmin") == cached_entry.get("dmax"):
                    if not do_cws or (cached_entry.get("cws_X") and cached_entry.get("cws_Z")):
                        if debug & 4:
                            print("[dist_m4ri] Cache hit for CSS distance (exact distance known)!")
                        dx_res = cached_entry.get("dX", format_bounds_list(cached_entry.get("dmin_X", 0), cached_entry.get("dmax_X", 0), cached_entry.get("rw_steps_X", 0)))
                        dz_res = cached_entry.get("dZ", format_bounds_list(cached_entry.get("dmin_Z", 0), cached_entry.get("dmax_Z", 0), cached_entry.get("rw_steps_Z", 0)))
                        return (
                            cached_entry["dist"], dx_res, dz_res,
                            cached_entry.get("cws_X", []), cached_entry.get("cws_Z", [])
                        ) if do_cws else (
                            cached_entry["dist"], dx_res, dz_res
                        )
                # Seed bounds from cache
                if eff_dmax == 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = cached_entry["dmax"]
                elif eff_dmax > 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = min(eff_dmax, cached_entry["dmax"])
                if eff_dmin <= 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = cached_entry["dmin"]
                elif eff_dmin > 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = max(eff_dmin, cached_entry["dmin"])
        except Exception:
            code_key = None
            cached_entry = None

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
        dmin_z, dmax_z, rw_steps_z = 0, 0, 0
        dmin_x, dmax_x, rw_steps_x = 0, 0, 0
        cws_Z, cws_X = [], []

        # Z-distance: Hx as finH, Hz as finG (or Lz as finL)
        if can_compute_Z:
            dmin_z, dmax_z, rw_steps_z = run_dist_m4ri(
                dist_m4ri_path=dist_m4ri,
                method=method,
                finH=file_Hx,
                finG=file_Hz if file_Lz is None else None,
                finL=file_Lz,
                dmin=eff_dmin,
                dmax=eff_dmax,
                wmin=wmin,
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
            if do_cws and outZ and os.path.exists(outZ):
                cws_Z = read_sparse_vectors(outZ)
                cws_Z.sort(key=len)

        # X-distance: Hz as finH, Hx as finG (or Lx as finL)
        if can_compute_X:
            dmin_x, dmax_x, rw_steps_x = run_dist_m4ri(
                dist_m4ri_path=dist_m4ri,
                method=method,
                finH=file_Hz,
                finG=file_Hx if file_Lx is None else None,
                finL=file_Lx,
                dmin=eff_dmin,
                dmax=eff_dmax,
                wmin=wmin,
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
            if do_cws and outX and os.path.exists(outX):
                cws_X = read_sparse_vectors(outX)
                cws_X.sort(key=len)

        dX_info = format_bounds_list(dmin_x, dmax_x, rw_steps_x) if can_compute_X else None
        dZ_info = format_bounds_list(dmin_z, dmax_z, rw_steps_z) if can_compute_Z else None

        if dist_X is not None and dist_Z is not None:
            dist = min(dist_Z, dist_X) if (dist_Z > 0 and dist_X > 0) else max(dist_Z, dist_X)
        elif dist_X is not None:
            dist = dist_X
        else:
            dist = dist_Z

        res_tuple = (dist, dX_info, dZ_info, cws_X, cws_Z) if do_cws else (dist, dX_info, dZ_info)
        if _use_distance_cache and code_key is not None:
            prev_steps = cached_entry.get("rw_steps", 0) if cached_entry else 0
            prev_dmax = cached_entry.get("dmax", 0) if cached_entry else 0
            prev_dmin = cached_entry.get("dmin", 0) if cached_entry else 0

            run_steps = (rw_steps_z if can_compute_Z else 0) + (rw_steps_x if can_compute_X else 0)
            total_rw_steps = prev_steps + run_steps

            curr_dmax = dist if dist > 0 else 0
            best_dmax = min(prev_dmax, curr_dmax) if (prev_dmax > 0 and curr_dmax > 0) else (curr_dmax if curr_dmax > 0 else prev_dmax)

            curr_dmin = 0
            if can_compute_Z and can_compute_X:
                curr_dmin = min(dmin_z, dmin_x)
            elif can_compute_Z:
                curr_dmin = dmin_z
            elif can_compute_X:
                curr_dmin = dmin_x
            best_dmin = max(prev_dmin, curr_dmin)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "dmin_X": dmin_x,
                "dmax_X": dmax_x,
                "rw_steps_X": rw_steps_x,
                "dmin_Z": dmin_z,
                "dmax_Z": dmax_z,
                "rw_steps_Z": rw_steps_z,
                "dX": dX_info,
                "dZ": dZ_info,
                "cws_X": cws_X,
                "cws_Z": cws_Z
            }
            if eff_cache_file:
                save_distance_cache(eff_cache_file)
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
    wmin: int = 1,
    wmax: int = 0,
    dW: int = -1,
    maxC: int = 0,
    pmin: float = 0.0,
    noscan: int = 0,
    do_cws: bool = False,
    cache_file: Optional[Union[str, Path]] = None,
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
        wmin: Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1).
        wmax: Maximum weight to search in CC.
        dW: Extra weight window above dmin to collect codewords.
        maxC: Maximum number of codewords to collect.
        pmin: Probability cutoff for error mechanisms in DEM.
        noscan: Skip CC scan loop if 1.
        do_cws: Whether to return extracted error mechanisms / codewords.
        cache_file: Optional JSON file path for persistent distance caching.
        solver: "dist_m4ri" or "codedistance".
        codedistance_method: Method if using codedistance library.
        codedistance_params: Extra parameters for codedistance library.
        seed: Random seed.
        debug: Debug level flags.

    Returns:
        tuple (dist, d_info, cws) if do_cws else (dist, d_info)
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
        if do_cws:
            raise ValueError("Codeword extraction (do_cws=True) is not supported with codedistance solver; use solver='dist_m4ri'.")
        codedistance = _get_codedistance()
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
        d_info = format_bounds_list(d, d, 0) if d > 0 else [0, 0, 0]
        return d, d_info

    # Solver is native multithreaded dist_m4ri
    global _distance_cache, _use_distance_cache, _distance_cache_file
    eff_cache_file = str(Path(cache_file).resolve()) if cache_file is not None else _distance_cache_file
    code_key = None
    cached_entry = None
    if _use_distance_cache:
        if eff_cache_file:
            load_distance_cache(eff_cache_file)
        try:
            dem_obj = dem if dem is not None else circuit
            dem_state = get_sparse_array_state(dem_obj)
            code_key = f"dem:{dem_state}" if pmin <= 0.0 else f"dem:{dem_state}:pmin={pmin}"
            cached_entry = _distance_cache.get(code_key)
            if cached_entry is not None:
                # If exact distance is already proven and not asking for more codewords
                if cached_entry.get("dmin", 0) > 0 and cached_entry.get("dmin") == cached_entry.get("dmax"):
                    if not do_cws or (cached_entry.get("cws") and len(cached_entry["cws"]) > 0):
                        if debug & 4:
                            print("[dist_m4ri] Cache hit for DEM distance (exact distance known)!")
                        d_info = format_bounds_list(cached_entry.get("dmin", 0), cached_entry.get("dmax", 0), cached_entry.get("rw_steps", 0))
                        return (
                            cached_entry["dist"], d_info, cached_entry.get("cws", [])
                        ) if do_cws else (
                            cached_entry["dist"], d_info
                        )
                # Seed bounds from cache
                if eff_dmax == 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = cached_entry["dmax"]
                elif eff_dmax > 0 and cached_entry.get("dmax", 0) > 0:
                    eff_dmax = min(eff_dmax, cached_entry["dmax"])
                if eff_dmin <= 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = cached_entry["dmin"]
                elif eff_dmin > 1 and cached_entry.get("dmin", 0) > 1:
                    eff_dmin = max(eff_dmin, cached_entry["dmin"])
        except Exception:
            code_key = None
            cached_entry = None

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

        dmin_res, dmax_res, rw_steps = run_dist_m4ri(
            dist_m4ri_path=dist_m4ri,
            method=method,
            fdem=file_dem,
            dmin=eff_dmin,
            dmax=eff_dmax,
            wmin=wmin,
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

        dist = dmin_res if (dmin_res == dmax_res or dmax_res == 0) else dmax_res
        cws = []
        if do_cws and outC and os.path.exists(outC):
            cws = read_sparse_vectors(outC)
            cws.sort(key=len)

        d_info = format_bounds_list(dmin_res, dmax_res, rw_steps)

        if _use_distance_cache and code_key is not None:
            prev_steps = cached_entry.get("rw_steps", 0) if cached_entry else 0
            prev_dmax = cached_entry.get("dmax", 0) if cached_entry else 0
            prev_dmin = cached_entry.get("dmin", 0) if cached_entry else 0
            prev_cws = list(cached_entry.get("cws", [])) if cached_entry else []

            total_rw_steps = prev_steps + rw_steps
            best_dmax = min(prev_dmax, dmax_res) if (prev_dmax > 0 and dmax_res > 0) else (dmax_res if dmax_res > 0 else prev_dmax)
            best_dmin = max(prev_dmin, dmin_res)

            combined_cws = prev_cws
            if cws:
                existing_set = {tuple(cw) for cw in combined_cws}
                for cw in cws:
                    if tuple(cw) not in existing_set:
                        combined_cws.append(cw)
                        existing_set.add(tuple(cw))
                combined_cws.sort(key=len)

            d_info = format_bounds_list(best_dmin, best_dmax, total_rw_steps)

            _distance_cache[code_key] = {
                "dist": dist,
                "dmin": best_dmin,
                "dmax": best_dmax,
                "rw_steps": total_rw_steps,
                "d_info": d_info,
                "cws": combined_cws
            }
            if eff_cache_file:
                save_distance_cache(eff_cache_file)

        if do_cws:
            return dist, d_info, cws
        return dist, d_info

    finally:
        for f in temp_files:
            if f and os.path.exists(f):
                try: os.remove(f)
                except OSError: pass


def _write_nzlist_file(filepath: str, cws: List[List[int]]) -> None:
    """Writes codewords to a text file in NZLIST format (1-based indices)."""
    with open(filepath, "w") as f:
        f.write("%% NZLIST\n")
        f.write(f"% {len(cws)} codewords\n")
        for cw in cws:
            f.write(f"{len(cw)} " + " ".join(str(idx + 1) for idx in cw) + "\n")


def parse_cli_args(argv: List[str]) -> Dict[str, Any]:
    """Parses CLI arguments supporting both key=value pairs and standard --flags."""
    args: Dict[str, Any] = {
        "method": 3,
        "finH": None,
        "finG": None,
        "finL": None,
        "fdem": None,
        "Hx": None,
        "Hz": None,
        "Lx": None,
        "Lz": None,
        "finC": None,
        "outC": None,
        "dmin": 0,
        "dmax": 0,
        "wmin": 1,
        "wmax": 0,
        "dexp": 0,
        "steps": None,
        "threads": None,
        "timeout": 60.0,
        "dW": -1,
        "maxC": 0,
        "pmin": 0.0,
        "noscan": 0,
        "classical": -1,
        "seed": 0,
        "debug": 0,
        "solver": "dist_m4ri",
        "cache_file": "tmp_dist_cache.json",
        "use_cache": True,
        "do_cws": False,
    }

    i = 0
    while i < len(argv):
        arg = argv[i]
        if not arg:
            i += 1
            continue

        if arg in ("-h", "--help", "help"):
            args["help"] = True
            i += 1
            continue

        if arg in ("--no-cache", "-no-cache", "nocache", "--nocache"):
            args["use_cache"] = False
            args["cache_file"] = None
            i += 1
            continue

        if arg in ("--cws", "-cws", "cws", "do_cws=1", "--do_cws"):
            args["do_cws"] = True
            i += 1
            continue

        key = None
        val = None
        if "=" in arg:
            key, val = arg.split("=", 1)
            if key.startswith("--"):
                key = key[2:]
            elif key.startswith("-"):
                key = key[1:]
        elif arg.startswith("--") or arg.startswith("-"):
            key = arg.lstrip("-")
            if i + 1 < len(argv) and not argv[i + 1].startswith("-") and "=" not in argv[i + 1]:
                val = argv[i + 1]
                i += 1
            else:
                val = "1"
        else:
            if os.path.exists(arg):
                if arg.endswith(".dem"):
                    args["fdem"] = arg
                elif arg.endswith(".mmx") or arg.endswith(".mtx"):
                    if args["finH"] is None:
                        args["finH"] = arg
                    elif args["finL"] is None:
                        args["finL"] = arg
            i += 1
            continue

        if key:
            key_lower = key.lower()
            if key_lower in ("cache", "cache_file", "cachefile"):
                if val.lower() in ("0", "none", "false", "off", "no"):
                    args["use_cache"] = False
                    args["cache_file"] = None
                else:
                    args["use_cache"] = True
                    args["cache_file"] = val
            elif key_lower in ("fdem", "dem"):
                args["fdem"] = val
            elif key_lower == "finh":
                args["finH"] = val
            elif key_lower == "fing":
                args["finG"] = val
            elif key_lower == "finl":
                args["finL"] = val
            elif key_lower in ("hx", "finhx"):
                args["Hx"] = val
            elif key_lower in ("hz", "finhz"):
                args["Hz"] = val
            elif key_lower in ("lx", "finlx"):
                args["Lx"] = val
            elif key_lower in ("lz", "finlz"):
                args["Lz"] = val
            elif key_lower == "finc":
                args["finC"] = val
            elif key_lower == "outc":
                args["outC"] = val
                args["do_cws"] = True
            elif key_lower in ("method", "m"):
                args["method"] = int(val)
            elif key_lower in ("dmin", "d_min"):
                args["dmin"] = int(val)
            elif key_lower in ("dmax", "d_max"):
                args["dmax"] = int(val)
            elif key_lower in ("wmin", "w_min"):
                args["wmin"] = int(val)
            elif key_lower in ("wmax", "w_max"):
                args["wmax"] = int(val)
            elif key_lower in ("dexp", "dest", "d_exp"):
                args["dexp"] = int(val)
            elif key_lower in ("steps", "num_steps", "nsteps"):
                args["steps"] = int(val)
            elif key_lower in ("threads", "num_threads", "t"):
                args["threads"] = int(val)
            elif key_lower in ("timeout", "time"):
                args["timeout"] = float(val)
            elif key_lower == "dw":
                args["dW"] = int(val)
            elif key_lower in ("maxc", "max_c"):
                args["maxC"] = int(val)
            elif key_lower == "pmin":
                args["pmin"] = float(val)
            elif key_lower == "noscan":
                args["noscan"] = int(val)
            elif key_lower == "classical":
                args["classical"] = int(val)
            elif key_lower == "seed":
                args["seed"] = int(val)
            elif key_lower in ("debug", "dbg"):
                args["debug"] = int(val)
            elif key_lower == "solver":
                args["solver"] = val
            elif key_lower == "cws":
                args["do_cws"] = bool(int(val)) if val.isdigit() else (val.lower() not in ("0", "false", "no"))

        i += 1

    return args


def print_cli_help() -> None:
    help_text = """dist_m4ri.py: Multithreaded distance calculator Python CLI

Usage: dist_m4ri.py [key=val | --flag val ...]

Options:
  fdem=FILE             Detector Error Model input file (.dem)
  finH=FILE             Parity check matrix input file (.mmx / .mtx)
  finG=FILE             Generator matrix input file (.mmx / .mtx)
  finL=FILE             Logical operator matrix input file (.mmx / .mtx)
  Hx=FILE, Hz=FILE      CSS check matrices (alternative to finH/finL)
  Lx=FILE, Lz=FILE      CSS logical operators (optional)
  method=N              1=RW, 2=CC, 3=Bracketing (default: 3)
  dmin=N                Certified lower bound, inclusive (default: 0)
  dmax=N                Known upper bound, inclusive (default: 0)
  wmin=N                Minimum distance of interest (terminate early if cw of weight <= wmin is found in RW or CC, default: 1)
  wmax=N                Maximum weight to search in CC
  dexp=N                Expected distance estimate
  steps=N               Maximum RW steps (default: 1000 in method 3)
  threads=N             Worker threads (default: hardware concurrency)
  timeout=SEC           Execution timeout in seconds (default: 60.0)
  dW=N                  Extra weight window above dmin to collect codewords
  maxC=N                Maximum number of codewords to collect
  outC=FILE             File to output non-trivial codewords (NZLIST format)
  pmin=PROB             Probability threshold for DEM errors
  noscan=1              Skip CC scan loop
  classical=1           Force classical mode (0 for CSS)
  seed=N                Random seed
  debug=N               Debug bitmask (e.g. 1, 2, 4)
  solver=NAME           'dist_m4ri' (default) or 'codedistance'
  cache=FILE            Persistent JSON cache file (default: tmp_dist_cache.json)
  --no-cache / nocache  Disable persistent JSON caching
  --cws                 Collect and output non-trivial codewords
"""
    print(help_text)


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    args = parse_cli_args(argv)

    if args.get("help") or (not args.get("fdem") and not args.get("finH") and not args.get("Hx")):
        print_cli_help()
        return 0

    cache_file = args["cache_file"] if args["use_cache"] else None
    if not args["use_cache"]:
        disable_distance_cache()

    try:
        if args["fdem"]:
            res = compute_dem_distance(
                dem=args["fdem"],
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                dW=args["dW"],
                maxC=args["maxC"],
                pmin=args["pmin"],
                noscan=args["noscan"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, d_info, cws = res
                print(format_bounds_str(d_info))
                if args["outC"]:
                    _write_nzlist_file(args["outC"], cws)
            else:
                dist, d_info = res
                print(format_bounds_str(d_info))
            return 0

        if args["Hx"] is not None or (args["classical"] == 0 and args["finH"] is not None):
            hx = args["Hx"] if args["Hx"] is not None else args["finH"]
            hz = args["Hz"] if args["Hz"] is not None else (args["finG"] if args["finG"] is not None else args["finL"])
            lx = args["Lx"]
            lz = args["Lz"]
            res = compute_css_distance(
                Hx=hx,
                Hz=hz,
                Lx=lx,
                Lz=lz,
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                dW=args["dW"],
                maxC=args["maxC"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, dx_info, dz_info, cws_x, cws_z = res
                exact_tag = " (exact)" if (dx_info and dx_info[0] > 0 and dx_info[0] == dx_info[1] and dz_info and dz_info[0] > 0 and dz_info[0] == dz_info[1]) else ""
                print(f"{dist} dX: {dx_info} dZ: {dz_info}{exact_tag}")
                if args["outC"]:
                    _write_nzlist_file(args["outC"], (cws_x or []) + (cws_z or []))
            else:
                dist, dx_info, dz_info = res
                exact_tag = " (exact)" if (dx_info and dx_info[0] > 0 and dx_info[0] == dx_info[1] and dz_info and dz_info[0] > 0 and dz_info[0] == dz_info[1]) else ""
                print(f"{dist} dX: {dx_info} dZ: {dz_info}{exact_tag}")
            return 0

        if args["finH"]:
            res = compute_classical_distance(
                H=args["finH"],
                method=args["method"],
                threads=args["threads"],
                timeout=args["timeout"],
                num_steps=args["steps"],
                d_exp=args["dexp"],
                dmin=args["dmin"],
                dmax=args["dmax"],
                wmin=args["wmin"],
                wmax=args["wmax"],
                dW=args["dW"],
                maxC=args["maxC"],
                do_cws=args["do_cws"] or (args["outC"] is not None),
                return_info=True,
                cache_file=cache_file,
                solver=args["solver"],
                seed=args["seed"],
                debug=args["debug"]
            )
            if args["do_cws"] or (args["outC"] is not None):
                dist, d_info, cws = res
                print(format_bounds_str(d_info))
                if args["outC"]:
                    _write_nzlist_file(args["outC"], cws)
            else:
                dist, d_info = res
                print(format_bounds_str(d_info))
            return 0
    except ValueError as e:
        sys.stderr.write(f"Error: {e}\n")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
