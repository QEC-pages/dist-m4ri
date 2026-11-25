"""Stub adapter for `tropidist` (implements Adapter interface with supports() == False).

Per REVIEW.md Section 3.3, this stub assumes no CLI or Python API for `tropidist`; implement
`prepare()`, `command()`, and `parse()` once `tropidist` is integrated and configured.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from . import Adapter


class TropidistAdapter(Adapter):
    """Stub adapter for `tropidist` (supports() is False until configured)."""

    name = "tropidist"
    interface = "cli"

    def __init__(self, config: str = "default", configured: bool = False) -> None:
        self.config = config
        self.configured = bool(configured)

    def describe(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "interface": self.interface,
            "config": self.config,
            "version": None,
            "commit": None,
            "dirty": None,
            "sha256": None,
            "deps": {},
            "features": {"stub": True, "configured": self.configured},
        }

    def supports(self, inst: Dict[str, Any]) -> bool:
        return False

    def prepare(self, inst: Dict[str, Any]) -> Dict[str, Any]:
        raise NotImplementedError("tropidist adapter is a stub; configure and implement after integration.")

    def command(
        self,
        inst: Dict[str, Any],
        threads: int,
        budget_s: float,
        seed: int,
        cfg: Dict[str, Any],
    ) -> Dict[str, Any]:
        raise NotImplementedError("tropidist adapter is a stub; configure and implement after integration.")

    def parse(
        self,
        rc: Optional[int],
        stdout: str,
        stderr: str,
        *,
        killed: bool = False,
        prep: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        raise NotImplementedError("tropidist adapter is a stub; configure and implement after integration.")
