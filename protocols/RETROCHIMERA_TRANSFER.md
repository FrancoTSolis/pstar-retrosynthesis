# Frozen RetroChimera transfer protocol

Written after a one-molecule adapter smoke test and before a RetroChimera
route search on a benchmark target.

## Models

Evaluate the public Pistachio and USPTO-FULL checkpoints with RetroChimera
1.3.0. Follow the model card by requesting top-10 candidates. Consensus mode,
which retains reactions proposed by both submodels, is the primary Pistachio
configuration. Default mode and top-50 are ablations.

## Planner

Transfer the ten frozen P* policies without changing any search parameter.
Only the single-step candidate generator is replaced. Keep the URSA stock,
N=500, depth 30, policy order, route deduplication, ten-route cap, and official
URSA/ChemCensor reranking.

## Benchmarks

Run in this order:

1. URSA-drugs&clinicals-2026
2. URSA-expert-2026
3. USPTO-190

Report P* Solv-0/1/2 and CC*, the included single policy, the paired portfolio
gain, wall time, and single-step calls. Published RetroChimera-MCTS values are
external comparisons rather than paired tests.
