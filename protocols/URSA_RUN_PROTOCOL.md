# Frozen URSA evaluation protocol

Written before C* or M* was run with the URSA stock or on the
URSA-drugs&clinicals-2026 targets.

## Versions and assets

- URSA tag `v1.1.0`, commit `b52a657`
- ChemCensor `v1.2.0`
- `ChemCensor-DB-U2-1.0.0.sqlite`
- `URSA_BBs_v1_0_0.csv`, Hugging Face revision
  `9ac1948b272fafbaaae99916a02eb5298907138b`
- `URSA-drugs-clinicals-2026.csv`, revision
  `616394dac893e1cf88509aaf77fadb2a813b4226`
- canonical-SMILES stock matching
- official URSA subtree collapsing and best-variant selection

## Frozen planners

- C*: full top-50 DESP view, building-block rule `2:1:0.5`,
  learned value alpha 16, width all
- M*: seeded random 20 of the DESP top 50, prior hidden,
  the same rule, learned value alpha 8, width all

The failed validity-retuning experiment does not alter either method.

## Runs

1. USPTO-190: C* and M*, N=500
2. URSA-drugs&clinicals-2026: C*, N=500

Search uses the official URSA stock and maximum route depth 30. Report
Solv-0/1/2 and CC* from official URSA, wall time, and hardware. No parameter
may change after either benchmark is observed.
