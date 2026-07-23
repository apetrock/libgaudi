"""Bracket and Newton search along scalar energy E(t)."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, List, Tuple

import numpy as np

EvalAll = Callable[[float], Tuple[float, float, float]]


@dataclass
class NewtonResult:
    t: float = 0.0
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    energy: float = 0.0
    energy_prime: float = 0.0
    residual: float = float("inf")
    iterations: int = 0
    converged: bool = False
    trace: List[Tuple[float, float]] = field(default_factory=list)


def bracket_min_energy(
    eval_all: EvalAll,
    t_min: float,
    t_max: float,
    n_samples: int = 200,
) -> Tuple[float, float]:
    best_t = t_min
    best_e = float("inf")
    for i in range(n_samples + 1):
        t = t_min + (t_max - t_min) * i / n_samples
        try:
            e, _, _ = eval_all(t)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            continue
        if e < best_e:
            best_e = e
            best_t = t
    return best_t, best_e


def newton_step(t: float, e_prime: float, e_pp: float, min_denom: float = 1e-14) -> float:
    """Undamped Newton step (used in unit tests)."""
    if abs(e_pp) < min_denom:
        return t
    return t - e_prime / e_pp


def damped_newton_step(
    t: float,
    e_prime: float,
    e_pp: float,
    *,
    min_denom: float = 1e-14,
    lm_lambda: float = 1e-8,
    grad_scale: float = 0.1,
) -> float:
    """1D Levenberg–Marquardt-style step when E'' is ill-conditioned."""
    if abs(e_pp) >= min_denom:
        return t - e_prime / e_pp
    if abs(e_prime) < 1e-14:
        return t
    return t - grad_scale * e_prime


def _clamp_t(t: float, t_max: float) -> float:
    return max(-t_max, min(t_max, t))


def _backoff_from_singularity(eval_all: EvalAll, t: float, t_max: float) -> Tuple[float, float, bool]:
    """If eval fails near |g|=0, return slightly interior point on same ray."""
    for frac in (0.999, 0.995, 0.99, 0.98, 0.95):
        t_try = math.copysign(abs(t) * frac, t)
        if abs(t_try) > t_max:
            continue
        try:
            e, e_p, _ = eval_all(t_try)
            return t_try, e, True
        except (ValueError, ZeroDivisionError, FloatingPointError):
            continue
    return t, float("inf"), False


def newton_min_energy(
    eval_all: EvalAll,
    t0: float,
    *,
    t_max: float,
    max_iters: int = 24,
    tol: float = 1e-6,
    min_denom: float = 1e-14,
    lm_lambda: float = 1e-8,
    stall_iters: int = 3,
    bracket_samples: int = 200,
) -> Tuple[float, float, bool]:
    """Hybrid Newton on E'(t)=0: damped NR, bracket fallback, backoff on eval failure."""
    t = _clamp_t(t0, t_max)
    e = float("inf")
    e_p = float("inf")
    stalls = 0

    for _ in range(max_iters):
        try:
            e, e_p, e_pp = eval_all(t)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            t_bo, e_bo, ok = _backoff_from_singularity(eval_all, t, t_max)
            if ok:
                return t_bo, e_bo, True
            break

        if abs(e_p) < tol:
            return t, e, True
        if abs(abs(t) - t_max) < 1e-4 * t_max and abs(e_p) < tol * 1000.0:
            return t, e, True

        t_new = _clamp_t(
            damped_newton_step(t, e_p, e_pp, min_denom=min_denom, lm_lambda=lm_lambda),
            t_max,
        )
        if abs(t_new - t) < tol * max(1.0, abs(t)):
            stalls += 1
        else:
            stalls = 0

        if stalls >= stall_iters:
            break
        t = t_new

    try:
        e, e_p, _ = eval_all(t)
    except (ValueError, ZeroDivisionError, FloatingPointError):
        t_bo, e_bo, ok = _backoff_from_singularity(eval_all, t, t_max)
        if ok:
            return t_bo, e_bo, True
        t_br, e_br = bracket_min_energy(eval_all, -t_max, t_max, n_samples=bracket_samples)
        if math.isfinite(e_br):
            return t_br, e_br, abs(e_p) < tol * 1000.0 if math.isfinite(e_p) else False
        return t, float("inf"), False

    if abs(e_p) < tol:
        return t, e, True
    if abs(abs(t) - t_max) < 0.05 * t_max:
        return t, e, True

    t_br, e_br = bracket_min_energy(eval_all, -t_max, t_max, n_samples=bracket_samples)
    if math.isfinite(e_br) and (not math.isfinite(e) or e_br <= e):
        return t_br, e_br, True
    return t, e, abs(e_p) < tol * 1000.0
