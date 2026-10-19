# P* Retrosynthesis

Reference implementation of the search methods from **Portfolio Search with
Precedent-Based Selection for Multi-Step Retrosynthesis**.

P* runs ten explicit policies in a shared AND–OR search harness. Each policy
uses the same single-step model and stock but applies a different combination
of reaction score and leaf value. The resulting routes are deduplicated and can
be sent to [URSA](https://github.com/insilicomedicine/URSA) for
ChemCensor-based precedent selection.

This repository is self-contained at the search layer. It includes:

- the AND–OR best-first engine;
- C*, M*, and the ten frozen P* policies;
- the exact fitted logistic leaf-value map used in the paper;
- stock, candidate-cache, route-validation, portfolio, and result interfaces;
- statistical utilities and unit tests;
- a toy data set that runs without model weights;
- adapters for cached predictions and public DESP assets.

Large model weights, reaction databases, and vendor catalogues are not
redistributed. The [DESP](https://github.com/coleygroup/desp),
[RetroChimera](https://github.com/microsoft/retrochimera), and
[URSA](https://github.com/insilicomedicine/URSA) projects publish the required
assets separately. Exact benchmark setup is documented in
[`docs/REPRODUCING.md`](docs/REPRODUCING.md).

## Installation

```bash
git clone git@github.com:FrancoTSolis/pstar-retrosynthesis.git
cd pstar-retrosynthesis
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Python 3.10 or newer is required. The base package uses RDKit, NumPy, PyYAML,
and pandas. Install the optional DESP dependencies with:

```bash
pip install -e ".[desp,dev]"
```

## Run the toy example

The toy example uses a small JSON candidate table and CSV stock. It exercises
the full search, all ten policies, route deduplication, and output validation.

```bash
pstar-portfolio \
  --config configs/pstar-toy.yaml \
  --targets examples/toy/targets.csv \
  --stock examples/toy/stock.csv \
  --candidates examples/toy/candidates.json \
  --output runs/toy

pstar-validate \
  --results runs/toy/portfolio.jsonl \
  --stock examples/toy/stock.csv \
  --candidates examples/toy/candidates.json
```

Run the test suite:

```bash
pytest
```

## Use cached model predictions

For a benchmark run, first generate or download one-step predictions. P* uses a
JSON mapping from canonical product SMILES to candidate records:

```json
{
  "CCO": [
    {
      "reactants": ["CC=O"],
      "score": 0.73,
      "template": "[C:1][O:2]>>[C:1]=[O:2]"
    }
  ]
}
```

Then run:

```bash
pstar-portfolio \
  --config configs/pstar.yaml \
  --targets data/targets.csv \
  --stock data/stock.csv \
  --candidates data/candidates.json \
  --output runs/benchmark
```

The search budget counts the first reveal of a distinct molecule within each
target. Reading a persistent prediction cache does not reduce this target-level
budget.

## Public DESP model

The optional DESP adapter reads the public `model_retro.pt` and
`idx2template_retro.json` files. Download the DESP assets from the source
project and run:

```bash
pstar-predict-desp \
  --model assets/model_retro.pt \
  --templates assets/idx2template_retro.json \
  --molecules data/molecules.csv \
  --output data/candidates.json \
  --top-n 50
```

The adapter applies templates with RDChiral, merges duplicate precursor sets,
and renormalizes their probabilities.

The exact logistic map is included. Exact leaf values also require the public
SCScore and DESP value outputs for each molecule. Supply these through the
`--descriptors` feature-table option as described in the reproduction guide.
The built-in fallbacks are intended only for the toy example.

## Official route-validity evaluation

This repository does not reimplement ChemCensor. Export the portfolio routes
and evaluate them with the official URSA release:

```bash
pstar-export-ursa \
  --results runs/benchmark/policy_*.jsonl \
  --output runs/benchmark/routes.json.gz
```

Use URSA v1.1.0, ChemCensor v1.2.0, the U2 database, and the pinned URSA
building-block catalogue to reproduce the paper's headline route-validity
results.

## Frozen policies

The definitions in `configs/pstar.yaml` are:

1. template prior + learned value, alpha 16
2. template prior + constant leaf value 2
3. no direct building-block reward + learned value, alpha 64
4. no direct building-block reward + learned value, alpha 8
5. template prior + DESP value
6. building-block rule + learned value, alpha 8
7. prior-gated building-block rule + learned value, alpha 16
8. C*: building-block rule + learned value, alpha 16
9. M*: matched random-20 view + building-block rule + learned value, alpha 8
10. template prior + SCScore value

The reaction rule is:

```text
2 * number_purchasable
- 1 * number_needing_synthesis
- 5 * largest_nonpurchasable_heavy_atom_ratio
```

## Repository layout

```text
configs/                 frozen YAML configurations
models/value_map.json    fitted nine-descriptor logistic map
protocols/               frozen experimental protocols
src/pstar/               installable Python package
tests/                   engine, policy, and portfolio tests
examples/toy/            model-free end-to-end example
```

## Reproducibility scope

The code and policy parameters are fully included. Full benchmark reproduction
also requires third-party model weights and data whose licenses do not permit
redistribution here. The paper Supporting Information records exact versions,
dataset revisions, statistical tests, and compute.

## License

MIT. See `LICENSE`.
