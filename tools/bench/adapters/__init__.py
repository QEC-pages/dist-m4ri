"""Engine adapters for the dist-m4ri certification benchmark (`tools/bench/`).

Each adapter implements the uniform interface of REVIEW.md Section 3.5:
  - name: str       ("dist_m4ri", "dist_m4ri.py", "qec_cert", "tropidist")
  - interface: str  ("binary", "python", "library", "cli")
  - config: str     ("default", "cc_fixed", "off", "auto", "cert", ...)
  - describe() -> dict
  - supports(inst) -> bool
  - prepare(inst) -> dict
  - command(inst, threads, budget_s, seed, cfg) -> dict
  - parse(rc, stdout, stderr, *, killed=False, prep=None) -> dict
"""
from __future__ import annotations

import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


class Adapter(ABC):
    """Abstract base class for benchmark engine adapters."""

    name: str
    interface: str
    config: str

    @abstractmethod
    def describe(self) -> Dict[str, Any]:
        """Return engine provenance: name, interface, config, version, commit, dirty, sha256, deps, features."""

    @abstractmethod
    def supports(self, inst: Dict[str, Any]) -> bool:
        """Return True iff this adapter configuration can run `inst`."""

    @abstractmethod
    def prepare(self, inst: Dict[str, Any]) -> Dict[str, Any]:
        """Convert or resolve input files (cached by sha256), returning at least {'prep_s': float}."""

    @abstractmethod
    def command(
        self,
        inst: Dict[str, Any],
        threads: int,
        budget_s: float,
        seed: int,
        cfg: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Build subprocess execution spec: {'argv': [...], 'env': {...}, 'cwd': str, 'argv_extra': [...]}."""

    @abstractmethod
    def parse(
        self,
        rc: Optional[int],
        stdout: str,
        stderr: str,
        *,
        killed: bool = False,
        prep: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Parse process output into the schema's partial `result` dict."""


def git_info(repo_dir: Path) -> Tuple[Optional[str], Optional[bool]]:
    """Return (short_commit, dirty_bool) for a git repository, or (None, None) if unavailable."""
    if not repo_dir.is_dir():
        return None, None
    try:
        head = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "--short=7", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if head.returncode != 0:
            return None, None
        commit = head.stdout.strip() or None
        st = subprocess.run(
            ["git", "-C", str(repo_dir), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        dirty = bool(st.stdout.strip()) if st.returncode == 0 else None
        return commit, dirty
    except Exception:
        return None, None


def build_adapter(
    spec: str,
    *,
    repo_root: Path,
    binary_path: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
    python_exe: Optional[str] = None,
    qec_cert_root: Optional[Path] = None,
    allow_fallback: bool = True,
) -> Adapter:
    """Instantiate an adapter from a spec string like 'dist_m4ri:default' or 'dist_m4ri.py:auto'."""
    from .dist_m4ri_bin import DistM4riBinAdapter
    from .dist_m4ri_py import DistM4riPyAdapter
    from .qec_cert import QecCertAdapter
    from .tropidist import TropidistAdapter

    if ":" in spec:
        name, config = spec.split(":", 1)
    else:
        name, config = spec, "default"

    if name == "dist_m4ri":
        return DistM4riBinAdapter(
            binary_path=binary_path,
            config=config,
            repo_root=repo_root,
            cache_dir=cache_dir,
            allow_fallback=allow_fallback,
        )
    if name == "dist_m4ri.py":
        return DistM4riPyAdapter(
            config=config,
            repo_root=repo_root,
            binary_path=binary_path,
            cache_dir=cache_dir,
            python_exe=python_exe,
            qec_cert_root=qec_cert_root,
        )
    if name == "qec_cert":
        return QecCertAdapter(config=config, qec_cert_root=qec_cert_root)
    if name == "tropidist":
        return TropidistAdapter(config=config)
    raise ValueError(f"Unknown adapter name {name!r} in spec {spec!r}")
