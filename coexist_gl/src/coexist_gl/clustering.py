from __future__ import annotations

import numpy as np


def _postorder_structure(tree, tip_order: list[str]):
    tip_index = {tip: i for i, tip in enumerate(tip_order)}
    nodes = list(tree.find_clades(order="postorder"))
    index = {id(node): i for i, node in enumerate(nodes)}
    children: list[tuple[int, ...]] = []
    leaf_positions: list[int] = []
    for node in nodes:
        children.append(tuple(index[id(child)] for child in node.clades))
        leaf_positions.append(tip_index[node.name] if not node.clades else -1)
    return children, leaf_positions


def fitch_steps(children: list[tuple[int, ...]], leaf_positions: list[int], labels: np.ndarray) -> int:
    masks = np.zeros(len(children), dtype=np.uint8)
    steps = 0
    for i, child_ids in enumerate(children):
        if not child_ids:
            masks[i] = 1 if labels[leaf_positions[i]] == 0 else 2
            continue
        mask = masks[child_ids[0]]
        for child in child_ids[1:]:
            intersection = mask & masks[child]
            if intersection:
                mask = intersection
            else:
                mask |= masks[child]
                steps += 1
        masks[i] = mask
    return steps


def clustering_test(tree, states: dict[str, int], permutations: int, seed: int) -> dict[str, object]:
    tip_order = [tip.name for tip in tree.get_terminals()]
    labels = np.array([states[tip] for tip in tip_order], dtype=np.uint8)
    children, leaf_positions = _postorder_structure(tree, tip_order)
    observed = fitch_steps(children, leaf_positions, labels)
    rng = np.random.default_rng(seed)
    at_least_as_clustered = 0
    for _ in range(permutations):
        permuted = rng.permutation(labels)
        at_least_as_clustered += fitch_steps(children, leaf_positions, permuted) <= observed
    p_value = (1 + at_least_as_clustered) / (permutations + 1)
    return {
        "fitch_steps_observed": observed,
        "permutation_reps": permutations,
        "p_cluster": p_value,
        "seed": seed,
        "status": "PHYLOGENETIC_CLUSTERING_SUPPORTED" if p_value < 0.05 else "PHYLOGENETIC_CLUSTERING_NOT_SUPPORTED",
    }

