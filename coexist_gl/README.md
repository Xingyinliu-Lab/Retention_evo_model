# coexist-gl v2 bundle

This is the current general host-level coexistence model. The formal data and model contract is `coexist_gl_family_v2` / `coexist_gl_v2_no_ale_unit_edge`. The Python package release string remains `0.1.0`; that release string does not denote an older scientific model.

For each host, the model observes `S=1 copy` or `C=2+ copies` on a rooted, binary, unit-edge species tree. It estimates the gain-to-resolution ratio `A=G/L` and tests the one-sided persistence parameter `r>=1`. `r=1` is the no-extra-persistence baseline; `r>1` increases the odds of `C->C` continuation.

## Included software

- `src/coexist_gl`: authoritative runtime implementation.
- `tests`: current tests.
- `scripts`: the final all-case/subset suite runner.
- `docs`: mathematical definition, input/output contract, installation and interpretation boundaries.

## Included cases

- `cases/actin_dalp_calp`: latest current584 actin case with the expanded `r` surface through 128.
- `cases/general_og_33`: the final 33 Archaea1338 de novo OGs. The earlier 31-case v2 outputs are preserved byte-for-byte; `OG0000323` and `OG0001582` were completed with the same current v2 runtime to close the final 33-case set.

AU and A/B labels are not inputs to coexist-gl. External reconciliation benchmarks are not packaged here.

## Run

Install locally with `python -m pip install -e .`, then:

```text
coexist-gl validate --config cases/actin_dalp_calp/input/family.yaml
coexist-gl run --config cases/actin_dalp_calp/input/family.yaml --output reproduced_actin_results
python scripts/run_example_suite.py --examples-root cases/general_og_33/cases --workers 8
```

Compare likelihoods only within the coexist-gl model and the same observation space.
