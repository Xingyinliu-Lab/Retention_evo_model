#!/usr/bin/env python3
"""Rebuild the frozen six-tree CTMC input from packaged source assets."""

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo


PACKAGE = Path(__file__).resolve().parents[1]
SOURCE_TABLES = PACKAGE / "input/source_tables"
SOURCE_TREES = PACKAGE / "input/source_trees"
PROFILE_META = SOURCE_TABLES / "frozen_authoritative_6769_seq_meta_with_aa_seq.tsv"
PLACEMENT = SOURCE_TABLES / "mapped_584.tsv"
EXPECTED = PACKAGE / "input/six_tree_tip_states.json"
PROFILE_MAP = {
    "canonical Actin-like": "cALP",
    "Crenactin-like / archaeal actin": "dALP",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def norm_tip(name):
    return name.split("_", 1)[1] if name.startswith(("GB_", "RS_")) else name


def build_states():
    counts = defaultdict(Counter)
    sequence_count = 0
    for row in read_rows(PROFILE_META):
        profile = PROFILE_MAP.get(row["actin_family_class"])
        if profile is None:
            continue
        counts[row["mag_id"]][profile] += 1
        sequence_count += 1
    if sequence_count != 2606 or len(counts) != 1662:
        raise ValueError(f"Frozen profile universe mismatch: sequences={sequence_count}, MAGs={len(counts)}")
    states = {}
    for mag in sorted(counts):
        current = counts[mag]
        states[mag] = (
            "coexist" if current["cALP"] and current["dALP"]
            else "cALP-only" if current["cALP"] else "dALP-only"
        )
    return states, sequence_count


def placement_concordant_states(states):
    labels = defaultdict(lambda: defaultdict(list))
    for row in read_rows(PLACEMENT):
        labels[row["mag_id"]][row["profile"]].append(row["branch_label"])
    retained = {}
    for mag, state in states.items():
        profiles = ("cALP", "dALP") if state == "coexist" else (
            ("cALP",) if state == "cALP-only" else ("dALP",)
        )
        concordant = True
        for profile in profiles:
            opposite = "dALP" if profile == "cALP" else "cALP"
            concordant &= profile + "_enriched" in labels[mag][profile]
            concordant &= opposite + "_enriched" not in labels[mag][profile]
        if concordant:
            retained[mag] = state
    return retained


def build_payload():
    states, sequence_count = build_states()
    retained = placement_concordant_states(states)
    trees = []
    audit_rows = []
    for path in sorted(SOURCE_TREES.glob("*.GTDB_rooted.treefile")):
        name = path.name.replace(".GTDB_rooted.treefile", "")
        tree = Phylo.read(str(path), "newick")
        tips = {norm_tip(tip.name or "") for tip in tree.get_terminals()}
        mapped = {mag: state for mag, state in retained.items() if mag in tips}
        if set(mapped) != tips:
            raise ValueError(f"Source tree is not closed against retained states: {name}")
        counts = Counter(mapped.values())
        trees.append({
            "name": name,
            "newick": path.read_text(encoding="utf-8"),
            "states": mapped,
            "root_semantics": "gtdb_projected",
        })
        audit_rows.append({
            "tree": name, "n_tips": len(tips),
            "dALP_only": counts["dALP-only"], "coexist": counts["coexist"],
            "cALP_only": counts["cALP-only"], "tree_sha256": sha256(path),
        })
    if len(trees) != 6:
        raise ValueError(f"Expected six source trees, observed {len(trees)}")
    payload = {
        "schema_version": "composition_corrected_five_model_input_v1",
        "state_encoding": {"dALP-only": 0, "coexist": 1, "cALP-only": 2},
        "tree_selection": "root_semantics == gtdb_projected",
        "trees": trees,
    }
    return payload, audit_rows, sequence_count, len(states), len(retained)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=PACKAGE / "audit/rebuilt_six_tree_tip_states.json")
    args = parser.parse_args()
    payload, audit_rows, sequence_count, mag_count, retained_count = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    expected_hash = sha256(EXPECTED)
    rebuilt_hash = sha256(args.output)
    audit = {
        "status": "PASS" if rebuilt_hash == expected_hash else "FAIL",
        "source_profile_sequences": sequence_count,
        "source_profile_MAGs": mag_count,
        "placement_concordant_MAGs": retained_count,
        "n_trees": len(audit_rows),
        "expected_sha256": expected_hash,
        "rebuilt_sha256": rebuilt_hash,
        "byte_identical": rebuilt_hash == expected_hash,
        "trees": audit_rows,
    }
    (PACKAGE / "audit/INPUT_REBUILD_AUDIT.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    if rebuilt_hash != expected_hash:
        raise SystemExit("Rebuilt input is not byte-identical to frozen input")
    print(json.dumps({"status": "PASS", "sha256": rebuilt_hash, "trees": len(audit_rows)}))


if __name__ == "__main__":
    main()
