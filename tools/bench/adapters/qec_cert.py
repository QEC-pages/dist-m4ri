"""Stub adapter for calling `qec-cert` directly (`qec_cert.pipeline`).

In v1 of the benchmark harness (dist-m4ri alone, ready for additional engines), `qec-cert` is
already benchmarked indirectly via `dist_m4ri.py:cert` and `dist_m4ri.py:auto`. Direct invocation
of `qec-cert` as a standalone engine is stubbed here and returns `supports() == False` until
explicitly configured and implemented.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
MANIFEST_DIR = REPO_ROOT / "benchmark" / "manifest"
if str(MANIFEST_DIR) not in sys.path:
    sys.path.insert(0, str(MANIFEST_DIR))
import make_inputs as mi  # noqa: E402

from . import Adapter, git_info  # noqa: E402
from .dist_m4ri_py import find_qec_cert_lib, find_qec_cert_root  # noqa: E402


def _qec_cert_version(qec_root: Optional[Path]) -> Optional[str]:
    if qec_root is None:
        return None
    for rel in ("qec_cert/__init__.py", "pyproject.toml"):
        p = qec_root / rel
        if not p.is_file():
            continue
        for line in p.read_text(encoding="utf-8").splitlines()[:60]:
            m = re.match(r'^(?:__version__|version)\s*=\s*["\']([^"\']+)["\']', line.strip())
            if m:
                return m.group(1)
    return None


class QecCertAdapter(Adapter):
    """Stub adapter for standalone `qec-cert` (supports() is False until configured)."""

    name = "qec_cert"
    interface = "library"

    def __init__(
        self,
        config: str = "default",
        qec_cert_root: Optional[Path] = None,
        configured: bool = False,
    ) -> None:
        self.config = config
        self.qec_cert_root = find_qec_cert_root(REPO_ROOT, explicit=qec_cert_root)
        self.configured = bool(configured)

    def describe(self) -> Dict[str, Any]:
        commit, dirty = git_info(self.qec_cert_root) if self.qec_cert_root else (None, None)
        lib_path = find_qec_cert_lib(self.qec_cert_root)
        return {
            "name": self.name,
            "interface": self.interface,
            "config": self.config,
            "version": _qec_cert_version(self.qec_cert_root),
            "commit": commit,
            "dirty": dirty,
            "sha256": mi.sha256_file(lib_path) if lib_path else None,
            "deps": {"libqeccert": str(lib_path) if lib_path else None},
            "features": {"stub": True, "configured": self.configured},
        }

    def supports(self, inst: Dict[str, Any]) -> bool:
        # Stubbed until direct qec-cert worker is wired up; use dist_m4ri.py:cert for now.
        return False

    def prepare(self, inst: Dict[str, Any]) -> Dict[str, Any]:
        raise NotImplementedError(
            "Direct qec_cert adapter is a stub in harness v1; use dist_m4ri.py:cert or implement direct worker."
        )

    def command(
        self,
        inst: Dict[str, Any],
        threads: int,
        budget_s: float,
        seed: int,
        cfg: Dict[str, Any],
    ) -> Dict[str, Any]:
        raise NotImplementedError("Direct qec_cert adapter is a stub in harness v1.")

    def parse(
        self,
        rc: Optional[int],
        stdout: str,
        stderr: str,
        *,
        killed: bool = False,
        prep: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        raise NotImplementedError("Direct qec_cert adapter is a stub in harness v1.")
