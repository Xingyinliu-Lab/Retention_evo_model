# Actin retention analysis reproducibility bundle

This directory contains the current, analysis-ready implementations and frozen case assets for two distinct analyses:

- `five_model/`: the composition-corrected three-state five-model comparison and ARD transition-contribution decomposition.
- `coexist_gl/`: the current host-level coexistence gain-loss program and its `coexist_gl_family_v2` cases.

It also includes two manuscript-level data assets:

- `6769_seq_meta_actin_sequences.fasta`: the 6,769 archaeal actin-superfamily protein sequences analysed in this study.
- `figure_source_data/`: one ZIP archive per manuscript figure, each containing its source data, plotting code and final assets.

The unpacked `source_data/` directory and the former monolithic figure-source-data ZIP are intentionally not tracked in this repository. Download the required figure-specific ZIP archive from `figure_source_data/`.

The two models answer different questions and their absolute likelihoods must not be compared. The five-model analysis uses `A-only / A+B coexist / B-only` or the corresponding actin profile states. `coexist-gl` ignores A/B labels and uses only family copy multiplicity per host (`S=1`, `C=2+`).

Included case sets are:

| Module | Case set |
|---|---|
| five-model | actin dALP/cALP, within-cALP A/B, MreB/ALP aggregate, 33 final general Archaea1338 OGs |
| coexist-gl | actin current584 dALP/cALP and the same 33 final general Archaea1338 OGs |

Historical automated1 within-cALP analyses, deprecated model variants, figures, logs, caches and AU intermediates are intentionally excluded.

Run `python tools/verify_package.py` from this directory to validate the bundled case counts, input closure and required result files.
