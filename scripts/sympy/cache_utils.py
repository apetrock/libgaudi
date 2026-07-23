"""Disk cache for expensive sympy derivations and lambdify."""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path
from typing import Any, Callable, TypeVar

import sympy as sp

CACHE_DIR = Path(__file__).resolve().parent / ".cache"

T = TypeVar("T")


def energy_cache_key() -> str:
    root = Path(__file__).resolve().parent
    h = hashlib.sha256()
    for name in ("shape_fcn.py", "compose.py"):
        h.update((root / name).read_bytes())
    h.update(sp.__version__.encode())
    return h.hexdigest()[:16]


def load_pickle(path: Path) -> Any:
    with path.open("rb") as f:
        return pickle.load(f)


def save_pickle(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as f:
        pickle.dump(value, f)


def disk_cache(path: Path, build: Callable[[], T]) -> T:
    if path.exists():
        try:
            return load_pickle(path)
        except (pickle.UnpicklingError, EOFError, AttributeError):
            pass
    value = build()
    save_pickle(path, value)
    return value
