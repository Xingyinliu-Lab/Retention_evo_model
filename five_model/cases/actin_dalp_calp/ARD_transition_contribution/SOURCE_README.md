# Composition-corrected ARD transition-route decomposition

## Status

Derived analysis of the frozen final five-model package. No evolutionary model was refitted and no frozen input was modified.

## Packaged inputs

- `../../input/six_tree_tip_states.json`
- `../../output/five_model_fits.tsv`

Only the 12 final `ARD_unrestricted` fits were used: six GTDB-projected species-tree reconstructions at unit-topology and source-branch-length scales, with the stationary distribution of the fitted Q matrix as the root prior.

## Definition

- `coexist-mediated`: dALP-only↔coexist and coexist↔cALP-only transitions.
- `direct endpoint`: dALP-only↔cALP-only transitions.

The values are fixed-Q posterior expected transition counts conditional on each frozen tree and its observed tip states. They are not literal duplication, transfer or loss event counts.

## Results and manuscript boundary

- Nine of 12 tree/scale histories were numerically stable under the predefined turnover/rate-boundary audit.
- Coexist-mediated routes contributed more than 50% of the posterior expected transition burden in 7/9 stable histories.
- The median direct-endpoint share among stable histories was 10.79%.
- Archaea1338 de novo was the stable exception: direct-endpoint shares were 53.52% with unit branches and 50.04% with source branches.
- Three source-branch histories were high-turnover or rate-boundary-saturated and are retained in the audit table but excluded from proportional interpretation.

Permitted conclusion: most numerically stable final ARD reconstructions reuse coexist-mediated adjacent routes, but this is not universal and the Archaea1338 de novo reconstruction supports an approximately balanced direct/coexist-mediated history under both branch scales.

Prohibited extrapolation: these route shares do not establish HGT, duplication/loss counts, cALP/dALP ancestry, absolute timing or a unique realized evolutionary history. The six trees are alternative reconstruction/sensitivity settings, not independent biological replicates.

## Outputs

- `figures/fig_ARD_transition_contribution_composition_corrected_20260901.pdf`
- `figures/fig_ARD_transition_contribution_composition_corrected_20260901.png`
- `source_data/ARD_transition_contribution_summary_composition_corrected.tsv`
- `source_data/ARD_transition_group_posterior_moments_composition_corrected.tsv`
- `source_data/ARD_individual_transition_posterior_moments_composition_corrected.tsv`
- `source_data/ARD_transition_contribution_composition_corrected_audit.json`
- `source_data/fig_ARD_transition_contribution_composition_corrected_20260901_source_data.tsv`
- `source_data/fig_ARD_transition_contribution_composition_corrected_20260901_excluded_histories.tsv`
- `code/decompose_composition_corrected_ARD.py`
- `code/plot_composition_corrected_ARD_contribution.R`
