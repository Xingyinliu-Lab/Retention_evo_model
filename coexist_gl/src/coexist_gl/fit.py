from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2

from .model import StateTreeLikelihood


LOG_A_BOUNDS = (-12.0, 12.0)
LOG_R_BOUNDS = (0.0, math.log(1_000_000.0))


@dataclass
class Point:
    r: float
    A: float
    log_likelihood: float


def _fit_A(model: StateTreeLikelihood, r: float) -> Point:
    result = minimize_scalar(
        lambda log_A: -model.evaluate(math.exp(log_A), r).log_likelihood,
        bounds=LOG_A_BOUNDS,
        method="bounded",
        options={"xatol": 1e-10},
    )
    A = math.exp(float(result.x))
    return Point(r, A, -float(result.fun))


def fit_null(model: StateTreeLikelihood) -> Point:
    return _fit_A(model, 1.0)


def fit_free(model: StateTreeLikelihood, null: Point, starts: int, seed: int) -> Point:
    rng = np.random.default_rng(seed)
    initial = [(math.log(null.A), x) for x in (0.0, 0.25, 1.0, 2.5)]
    while len(initial) < starts:
        initial.append((math.log(null.A) + rng.normal(0, 1), rng.uniform(0, 4)))
    best = Point(1.0, null.A, null.log_likelihood)
    for start in initial[:starts]:
        result = minimize(
            lambda x: -model.evaluate(math.exp(float(x[0])), math.exp(float(x[1]))).log_likelihood,
            np.asarray(start),
            method="L-BFGS-B",
            bounds=[LOG_A_BOUNDS, LOG_R_BOUNDS],
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 1000},
        )
        candidate = Point(math.exp(float(result.x[1])), math.exp(float(result.x[0])), -float(result.fun))
        if candidate.log_likelihood > best.log_likelihood:
            best = candidate
    return best


def profile_r(model: StateTreeLikelihood, free: Point) -> list[Point]:
    maximum = max(64.0, free.r * 16.0)
    grid = np.unique(np.concatenate(([1.0, free.r], np.exp(np.linspace(0.0, math.log(maximum), 81)))))
    return [_fit_A(model, float(r)) for r in grid]


def fit_model(model: StateTreeLikelihood, starts: int, seed: int) -> dict[str, object]:
    null = fit_null(model)
    free = fit_free(model, null, starts, seed)
    delta = max(0.0, free.log_likelihood - null.log_likelihood)
    statistic = 2.0 * delta
    p_value = 1.0 if statistic <= 0 else 0.5 * float(chi2.sf(statistic, 1))
    profile = profile_r(model, free)
    cutoff = free.log_likelihood - 0.5 * float(chi2.ppf(0.95, 1))
    accepted = [point.r for point in profile if point.log_likelihood >= cutoff]
    upper_open = bool(accepted and max(accepted) == max(point.r for point in profile))
    identifiable = not upper_open
    if not identifiable:
        status = "R_NOT_IDENTIFIABLE"
    elif free.r > 1.0 + 1e-7 and p_value < 0.05:
        status = "COEXIST_PERSISTENCE_SUPPORTED"
    else:
        status = "COEXIST_OBSERVED_NULL_BOUNDARY"
    return {
        "null": null,
        "free": free,
        "delta_log_likelihood": delta,
        "lrt_statistic": statistic,
        "p_persistence_one_sided": p_value,
        "r_profile_low": min(accepted) if accepted else None,
        "r_profile_high": None if upper_open or not accepted else max(accepted),
        "status": status,
        "profile": profile,
    }


def surface_A_values(free_A: float, configured: tuple[float, ...] | None) -> list[float]:
    if configured is not None:
        return sorted(set(configured))
    factors = [2 ** x for x in (-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3)]
    return sorted(set(free_A * factor for factor in factors))

