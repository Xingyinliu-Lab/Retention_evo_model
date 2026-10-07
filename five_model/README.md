# Composition-corrected five-model analysis

## Current model set

The current implementation compares five continuous-time three-state models on a rooted species tree:

1. `ARD_unrestricted`;
2. `adjacency_coexist_middle` (Sequential Retention);
3. `adjacency_B_middle`;
4. `adjacency_A_middle`;
5. `ancestral_coexistence_strict`.

For the reversible models, the root prior is the stationary distribution implied by the fitted rate matrix. The strict ancestral-coexistence model fixes the root state to coexistence. Fits are performed separately for unit topology and source branch lengths. AIC and delta AIC are calculated only within the same case and branch treatment.

The reference implementation is under `engine/composition_corrected/`. Case-specific wrappers are retained where the input representation differs.

The ARD transition-contribution code is retained separately under `engine/ard_transition_contribution/`. Its corresponding frozen outputs are included for every packaged case set; it is not merely a plotting step.

## Case directories

- `cases/actin_dalp_calp`: frozen six-tree dALP/cALP analysis.
- `cases/actin_within_calp`: current **gappyout584-backbone** within-cALP A/B analysis; no automated1 result is included.
- `cases/actin_mreb_alp`: MreB/Mbl-profile aggregate versus ALP aggregate analysis.
- `cases/general_og_33`: final 33 Retention-classified Archaea1338 de novo OG cases, each with frozen host states, induced GTDB r226 species tree, unit outputs and source-branch outputs.

The 33-case set is defined by the final Retention classification; AU is orthogonal evidence and is not a membership gate.

## General-OG execution

From `cases/general_og_33`, run the unit, ARD and source-branch analyses using the bundled scripts. The input contract for every case is `input/host_states.tsv`, `input/species_tree.nwk` and `input/INPUT_FREEZE.json`. See `run_all.ps1` or `run_all.sh`.
