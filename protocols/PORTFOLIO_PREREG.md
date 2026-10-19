# Frozen P* portfolio protocol

Written before any portfolio policy was run with the URSA stock or combined
into an official URSA best-of-10 evaluation.

## Motivation

On the original dev_a routes, the best single code policy reaches Solv-2 on
28 of 200 targets and the union reaches 55. Greedy set cover reaches all 55
with seven policies. This licensed a portfolio test but did not establish
performance under the smaller URSA stock.

## Frozen policies

1. template prior + learned value, alpha 16
2. template prior + constant leaf value 2
3. no direct building-block reward + learned value, alpha 64
4. no direct building-block reward + learned value, alpha 8
5. template prior + DESP value, scale 4
6. building-block rule + learned value, alpha 8
7. prior-gated building-block rule, threshold 0.02, + learned value, alpha 16
8. original C*: building-block rule + learned value, alpha 16
9. M*: matched-view seed 0 + building-block rule + learned value, alpha 8
10. template prior + SCScore value, scale 1

Each policy emits at most one route at N=500. Duplicate routes are removed.
The remaining routes, in the fixed order above, are passed to official URSA.

## Development gate

- Run all policies on the first 200 dev_a targets with the URSA stock.
- Compare P* with policy 8 on the same targets.
- Confirm the unchanged portfolio once on the first 200 dev_b targets.
- Proceed to final benchmarks only if Solv-2 improves by at least five
  percentage points on both sets.

No policy may be added, removed, reordered, or reweighted after an URSA-stock
development result is observed.
