# ARD transition-contribution implementations

The ARD decomposition is part of the final five-model workflow. It conditions on the fitted unrestricted ARD rate matrix and computes posterior expected directed-jump contributions on the same species tree and tip states.

Two groups are reported:

- `coexist_mediated`: the four endpoint-to-coexistence or coexistence-to-endpoint transitions;
- `direct_endpoint`: the two direct A-only to B-only transitions.

Included implementations:

- `decompose_six_tree_composition_corrected_ARD.py`: the final six-tree dALP/cALP implementation.
- `compute_ard_contribution_one_case.py`: the final one-case implementation used for within-cALP and general-OG cases. It accepts `--branch-scale unit_topology` (default) or `--branch-scale source_branch_lengths`, reads the matched ARD fit for that scale, and writes scale-specific contribution outputs.

Examples:

```bash
python compute_ard_contribution_one_case.py /path/to/case --branch-scale unit_topology
python compute_ard_contribution_one_case.py /path/to/case --branch-scale source_branch_lengths
```

MreB/ALP-specific state labels use the equivalent implementation retained under `cases/actin_mreb_alp/code/`.

These values are posterior expected transition contributions, not literal duplication, transfer or loss counts.
