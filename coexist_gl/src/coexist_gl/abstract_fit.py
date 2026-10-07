from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

from .abstract_model import AbstractTreeLikelihood


def fit_abstract(model: AbstractTreeLikelihood, starts: int, seed: int) -> dict[str, object]:
    rng = np.random.default_rng(seed)

    def null_objective(x):
        return -model.evaluate(math.exp(float(x[0])), 1.0, math.exp(float(x[1]))).log_likelihood

    null_starts = [(math.log(0.5), math.log(10.0)), (math.log(0.1), 0.0), (math.log(2.0), math.log(100.0))]
    null_results = [minimize(null_objective, x, method="L-BFGS-B", bounds=[(-12, 12), (-6, 14)]) for x in null_starts]
    null = min(null_results, key=lambda x: x.fun)

    def free_objective(x):
        return -model.evaluate(math.exp(float(x[0])), math.exp(float(x[1])), math.exp(float(x[2]))).log_likelihood

    free_starts = [(float(null.x[0]), theta, float(null.x[1])) for theta in (0.0, 0.5, 2.0, 4.0)]
    while len(free_starts) < starts:
        free_starts.append((float(null.x[0]) + rng.normal(), rng.uniform(0, 5), float(null.x[1]) + rng.normal()))
    free_results = [
        minimize(free_objective, x, method="L-BFGS-B", bounds=[(-12, 12), (0, math.log(1_000_000)), (-6, 14)])
        for x in free_starts[:starts]
    ]
    free = min(free_results, key=lambda x: x.fun)
    ll0, ll1 = -float(null.fun), -float(free.fun)
    statistic = max(0.0, 2.0 * (ll1 - ll0))
    return {
        "null_r": 1.0,
        "null_A": math.exp(float(null.x[0])),
        "null_kappa": math.exp(float(null.x[1])),
        "null_logL": ll0,
        "free_r": math.exp(float(free.x[1])),
        "free_A": math.exp(float(free.x[0])),
        "free_kappa": math.exp(float(free.x[2])),
        "free_logL": ll1,
        "delta_logL": ll1 - ll0,
        "LRT_stat": statistic,
        "p_persistence_one_sided": 1.0 if statistic <= 0 else 0.5 * float(chi2.sf(statistic, 1)),
        "interpretation_scope": "coarse_species_cut_sensitivity",
    }

