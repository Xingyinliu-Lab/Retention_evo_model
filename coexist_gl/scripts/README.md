# Final suite runner

The authoritative model runtime is `../src/coexist_gl/`. This directory retains
only `run_example_suite.py`, which can execute all final cases or selected cases.

Run all 33 final general-OG cases:

```bash
python scripts/run_example_suite.py --examples-root cases/general_og_33/cases --workers 8
```

Run one or more selected cases by adding repeated `--case CASE_ID` arguments.
No historical benchmark scripts or outputs are included in this module.
