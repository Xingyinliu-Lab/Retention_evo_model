from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class AbstractNode:
    name: str
    children: list["AbstractNode"] = field(default_factory=list)
    n_single: int = 0
    n_coexist: int = 0


@dataclass
class SpeciesCut:
    cut_id: str
    k: int
    root: AbstractNode
    members: list[dict[str, object]]


def _leaf_count(node) -> int:
    return len(node.get_terminals())


def discover_topology_cut(tree, states: dict[str, int], k: int) -> SpeciesCut:
    n_tips = len(tree.get_terminals())
    if k < 2 or k > n_tips:
        raise ValueError(f"topology_k must be between 2 and {n_tips}")
    frontier = [tree.root]
    while len(frontier) < k:
        candidates = [(-_leaf_count(node), i, node) for i, node in enumerate(frontier) if node.clades]
        if not candidates:
            raise ValueError(f"tree cannot be divided into {k} topology blocks")
        _, index, node = min(candidates)
        if len(frontier) - 1 + len(node.clades) > k:
            alternatives = [x for x in candidates if len(frontier) - 1 + len(x[2].clades) <= k]
            if not alternatives:
                raise ValueError(f"tree cannot reach exactly {k} blocks")
            _, index, node = min(alternatives)
        frontier[index:index + 1] = list(node.clades)

    chosen = {id(node): f"block_{i+1:04d}" for i, node in enumerate(frontier)}
    members: list[dict[str, object]] = []
    for node in frontier:
        block = chosen[id(node)]
        members.extend(
            {"species_cut_id": f"topology_k{k}", "species_block": block, "host_id": tip.name}
            for tip in sorted(node.get_terminals(), key=lambda x: x.name)
        )

    internal_counter = 0

    def contract(node) -> AbstractNode:
        nonlocal internal_counter
        if id(node) in chosen:
            tips = [tip.name for tip in node.get_terminals()]
            return AbstractNode(
                name=chosen[id(node)],
                n_single=sum(states[tip] == 0 for tip in tips),
                n_coexist=sum(states[tip] == 1 for tip in tips),
            )
        internal_counter += 1
        return AbstractNode(name=f"internal_{internal_counter:04d}", children=[contract(child) for child in node.clades])

    return SpeciesCut(f"topology_k{k}", k, contract(tree.root), members)

