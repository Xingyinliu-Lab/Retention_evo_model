from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


def _logadd(a: float, b: float) -> float:
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    high = max(a, b)
    return high + math.log1p(math.exp(min(a, b) - high))


def transition_matrix(A: float, r: float) -> np.ndarray:
    if A <= 0 or r < 1:
        raise ValueError("A must be >0 and r must be >=1")
    decay = np.exp(-(A + 2.0))
    pi_s = 2.0 / (A + 2.0)
    pi_c = A / (A + 2.0)
    base = np.array(
        [
            [pi_s + pi_c * decay, pi_c * (1.0 - decay)],
            [pi_s * (1.0 - decay), pi_c + pi_s * decay],
        ],
        dtype=float,
    )
    tilted = base.copy()
    denominator = base[1, 0] + r * base[1, 1]
    tilted[1, 1] = r * base[1, 1] / denominator
    tilted[1, 0] = 1.0 - tilted[1, 1]
    return tilted


def root_prior(A: float) -> np.ndarray:
    return np.array([2.0 / (A + 2.0), A / (A + 2.0)], dtype=float)


@dataclass
class LikelihoodResult:
    log_likelihood: float
    expected: dict[str, float] | None = None
    branch_rows: list[dict[str, object]] | None = None


class StateTreeLikelihood:
    def __init__(self, tree, tip_states: dict[str, int]):
        self.tree = tree
        self.tip_states = tip_states
        self.postorder = list(tree.find_clades(order="postorder"))
        self.preorder = list(tree.find_clades(order="preorder"))
        self.node_index = {id(node): i for i, node in enumerate(self.postorder)}
        self.children = [tuple(self.node_index[id(child)] for child in node.clades) for node in self.postorder]
        self.leaf_states = np.full(len(self.postorder), -1, dtype=np.int8)
        self.node_labels: dict[int, str] = {}
        internal_index = 0
        for node in self.preorder:
            if node.clades:
                internal_index += 1
                self.node_labels[id(node)] = f"N{internal_index:04d}"
            else:
                self.node_labels[id(node)] = str(node.name)
                self.leaf_states[self.node_index[id(node)]] = tip_states[str(node.name)]
        self.preorder_indices = [self.node_index[id(node)] for node in self.preorder]

    def evaluate(self, A: float, r: float, with_expectations: bool = False) -> LikelihoodResult:
        P = transition_matrix(A, r)
        logP = np.log(P)
        inside = np.zeros((len(self.postorder), 2), dtype=float)
        child_contribution: dict[tuple[int, int], tuple[float, float]] = {}

        for node_i, child_ids in enumerate(self.children):
            if not child_ids:
                state = self.leaf_states[node_i]
                inside[node_i, :] = (0.0, -math.inf) if state == 0 else (-math.inf, 0.0)
                continue
            total0 = 0.0
            total1 = 0.0
            for child_i in child_ids:
                contribution = (
                    _logadd(logP[0, 0] + inside[child_i, 0], logP[0, 1] + inside[child_i, 1]),
                    _logadd(logP[1, 0] + inside[child_i, 0], logP[1, 1] + inside[child_i, 1]),
                )
                child_contribution[(node_i, child_i)] = contribution
                total0 += contribution[0]
                total1 += contribution[1]
            inside[node_i, :] = (total0, total1)

        log_pi = np.log(root_prior(A))
        root_i = self.node_index[id(self.tree.root)]
        log_likelihood = _logadd(log_pi[0] + inside[root_i, 0], log_pi[1] + inside[root_i, 1])
        if not with_expectations:
            return LikelihoodResult(log_likelihood)

        outside = np.full((len(self.postorder), 2), -math.inf, dtype=float)
        outside[root_i, :] = log_pi
        expected = np.zeros((2, 2), dtype=float)
        branch_rows: list[dict[str, object]] = []
        for parent in self.preorder:
            parent_i = self.node_index[id(parent)]
            child_ids = self.children[parent_i]
            if not child_ids:
                continue
            for child in parent.clades:
                child_i = self.node_index[id(child)]
                sibling0 = 0.0
                sibling1 = 0.0
                for sibling_i in child_ids:
                    if sibling_i != child_i:
                        contribution = child_contribution[(parent_i, sibling_i)]
                        sibling0 += contribution[0]
                        sibling1 += contribution[1]
                outside[child_i, 0] = _logadd(
                    outside[parent_i, 0] + sibling0 + logP[0, 0],
                    outside[parent_i, 1] + sibling1 + logP[1, 0],
                )
                outside[child_i, 1] = _logadd(
                    outside[parent_i, 0] + sibling0 + logP[0, 1],
                    outside[parent_i, 1] + sibling1 + logP[1, 1],
                )
                sibling_sum = np.array([sibling0, sibling1])
                joint_log = (
                    outside[parent_i, :][:, None]
                    + sibling_sum[:, None]
                    + logP
                    + inside[child_i, :][None, :]
                    - log_likelihood
                )
                joint = np.exp(joint_log)
                joint /= joint.sum()
                expected += joint
                branch_rows.append(
                    {
                        "parent": self.node_labels[id(parent)],
                        "child": self.node_labels[id(child)],
                        "posterior_S_to_S": joint[0, 0],
                        "posterior_S_to_C": joint[0, 1],
                        "posterior_C_to_S": joint[1, 0],
                        "posterior_C_to_C": joint[1, 1],
                    }
                )

        named = {
            "expected_S_to_S": float(expected[0, 0]),
            "expected_S_to_C_gain": float(expected[0, 1]),
            "expected_C_to_S_resolution": float(expected[1, 0]),
            "expected_C_to_C_persistence": float(expected[1, 1]),
            "expected_changed_edges": float(expected[0, 1] + expected[1, 0]),
            "expected_total_edges": float(expected.sum()),
        }
        return LikelihoodResult(log_likelihood, named, branch_rows)
