# Final audit

- Formal scope: five composition-corrected three-state CTMC models, six species-tree reconstructions, unit topology.
- Fits: 30/30 present, finite and optimizer-successful.
- Frozen state counts: GTDB uses 4,208 exact-mapped MAGs; each A1338 tree uses 284 MAGs (64 M-only, 31 coexist, 189 ALP-only); each Asgard971 tree uses 851 MAGs (14, 182, 655).
- Root policy: stationary(Q) for ARD and all three reversible adjacency graphs; fixed coexist root only for strict ancestral coexistence.
- Best model: ARD for GTDB and all three A1338 trees; sequential retention for both Asgard971 trees.
- Paired against sequential retention across the six matched tree reconstructions: ARD median AIC difference = -18.755 (2/6 favor retention; exact two-sided P=0.15625); ALP-middle = 27.464 (3/6; P=0.4375); M-middle = 19.872 (5/6; P=0.15625); strict ancestral coexistence = 74.514 (6/6; P=0.03125).

Interpretation: strict ancestral coexistence is consistently disfavored. Sequential retention is strongly supported in the Asgard971 sensitivity layer but is not the universal best model across GTDB and Archaea1338, where unrestricted ARD fits better. This test does not infer evolutionary direction or D/T/L events.
