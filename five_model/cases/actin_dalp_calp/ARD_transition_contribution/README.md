# dALP/cALP ARD transition-contribution analysis

This directory contains the final composition-corrected ARD route decomposition dated 2026-09-01.

Inputs are read from the sibling frozen five-model directories:

- `../input/six_tree_tip_states.json`;
- `../output/five_model_fits.tsv`;
- `../code/fit_composition_corrected_five_models.py`.

Run:

```text
python code/decompose_composition_corrected_ARD.py
```

The analysis generates individual transition moments, grouped transition moments, the coexist-mediated/direct summary and a numerical audit under `output/`.

