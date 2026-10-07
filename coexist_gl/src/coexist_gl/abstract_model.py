from __future__ import annotations

import math

import numpy as np
from scipy.special import betaln, gammaln

from .model import LikelihoodResult, _logadd, root_prior, transition_matrix
from .species_cuts import AbstractNode


def _log_beta_binomial(n: int, m: int, p: float, kappa: float) -> float:
    alpha = p * kappa
    beta = (1.0 - p) * kappa
    return float(
        gammaln(n + 1)
        - gammaln(m + 1)
        - gammaln(n - m + 1)
        + betaln(m + alpha, n - m + beta)
        - betaln(alpha, beta)
    )


class AbstractTreeLikelihood:
    def __init__(self, root: AbstractNode):
        self.root = root
        self.postorder: list[AbstractNode] = []

        def visit(node: AbstractNode) -> None:
            for child in node.children:
                visit(child)
            self.postorder.append(node)

        visit(root)
        self.index = {id(node): i for i, node in enumerate(self.postorder)}

    def evaluate(self, A: float, r: float, kappa: float, with_expectations: bool = False) -> LikelihoodResult:
        if kappa <= 0:
            raise ValueError("kappa must be >0")
        P = transition_matrix(A, r)
        logP = np.log(P)
        inside = np.zeros((len(self.postorder), 2), dtype=float)
        for node in self.postorder:
            i = self.index[id(node)]
            if not node.children:
                n = node.n_single + node.n_coexist
                inside[i, 0] = _log_beta_binomial(n, node.n_coexist, P[0, 1], kappa)
                inside[i, 1] = _log_beta_binomial(n, node.n_coexist, P[1, 1], kappa)
                continue
            total = np.zeros(2, dtype=float)
            for child in node.children:
                child_i = self.index[id(child)]
                total[0] += _logadd(logP[0, 0] + inside[child_i, 0], logP[0, 1] + inside[child_i, 1])
                total[1] += _logadd(logP[1, 0] + inside[child_i, 0], logP[1, 1] + inside[child_i, 1])
            inside[i, :] = total
        pi = root_prior(A)
        root_i = self.index[id(self.root)]
        log_likelihood = _logadd(math.log(pi[0]) + inside[root_i, 0], math.log(pi[1]) + inside[root_i, 1])
        return LikelihoodResult(log_likelihood)

