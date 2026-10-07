from __future__ import annotations

import hashlib


def discover_gene_splits(gene_tree) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    all_tips = frozenset(tip.name for tip in gene_tree.get_terminals())
    seen: set[frozenset[str]] = set()
    rows: list[dict[str, object]] = []
    members: list[dict[str, object]] = []
    for node in gene_tree.find_clades(order="preorder"):
        if node is gene_tree.root or not node.clades:
            continue
        side = frozenset(tip.name for tip in node.get_terminals())
        other = all_tips - side
        if len(side) < 2 or len(other) < 2:
            continue
        canonical = side if (len(side), sorted(side)[0]) <= (len(other), sorted(other)[0]) else other
        if canonical in seen:
            continue
        seen.add(canonical)
        member_text = "\n".join(sorted(canonical)).encode()
        support = node.confidence if node.confidence is not None else "NA"
        split_id = f"split_{len(rows)+1:04d}"
        rows.append(
            {
                "split_id": split_id,
                "n_genes_side": len(canonical),
                "n_genes_other": len(all_tips) - len(canonical),
                "branch_support": support,
                "member_sha256": hashlib.sha256(member_text).hexdigest(),
                "first_gene": min(canonical),
            }
        )
        members.extend({"split_id": split_id, "gene_id": gene} for gene in sorted(canonical))
    return rows, members
