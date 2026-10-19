# Reproducing the paper experiments

## What this repository provides

The search engine, all reaction scorers, the exact fitted logistic value map,
the frozen portfolio, route validation, statistical utilities, protocol
documents, and aggregate result tables are versioned here.

## Third-party assets

Download these from their source projects:

1. **DESP**
   - `model_retro.pt`
   - `idx2template_retro.json`
   - `retro_value.pt`
   - the eMolecules building-block catalogue used by DESP
2. **RetroChimera 1.3.0**
   - Pistachio checkpoint
   - USPTO-FULL checkpoint
3. **URSA 1.1.0**
   - `URSA_BBs_v1_0_0.csv`
   - URSA benchmark target CSV files
4. **ChemCensor 1.2.0**
   - `ChemCensor-DB-U2-1.0.0.sqlite`

The paper used the dataset revisions recorded in
`protocols/URSA_RUN_PROTOCOL.md`.

## Exact value features

`models/value_map.json` is the exact fitted logistic map. Its nine inputs are:

1. heavy atoms divided by 10
2. total ring count
3. aromatic ring count
4. stereocenter count
5. heteroatom count divided by 5
6. rotatable bond count divided by 5
7. SCScore
8. synthetic-accessibility score divided by 3
9. non-negative DESP value divided by 3

For exact benchmark reproduction, supply a JSON feature table through
`--descriptors`, or connect the public SCScore and DESP value models to
`PaperFeatureProvider`. If no table is supplied, the CLI uses documented
fallbacks for the two external learned descriptors. Those fallbacks are for the
toy example and are not claimed to reproduce paper numbers.

## Route finding

Generate DESP candidates, then run C*, M*, or the full portfolio:

```bash
pstar-predict-desp \
  --model assets/model_retro.pt \
  --templates assets/idx2template_retro.json \
  --molecules data/molecules.csv \
  --output data/candidates.json \
  --top-n 50 \
  --device cuda:0

pstar-search \
  --config configs/cstar.yaml \
  --policy P08_Cstar \
  --targets data/targets.csv \
  --stock data/emolecules.csv \
  --candidates data/candidates.json \
  --descriptors data/descriptors.json \
  --output runs/cstar.jsonl
```

## Portfolio and URSA

```bash
pstar-portfolio \
  --config configs/pstar.yaml \
  --targets data/targets.csv \
  --stock data/ursa_stock.csv \
  --candidates data/candidates.json \
  --descriptors data/descriptors.json \
  --output runs/pstar

pstar-export-ursa \
  --results runs/pstar/P*.jsonl \
  --output runs/pstar/routes.json.gz
```

Run the exported routes with official URSA using top-k 10 and the pinned stock.
URSA performs the collapsed-route generation and final precedent selection.

## Verification

```bash
pstar-validate \
  --results runs/pstar/portfolio.jsonl \
  --stock data/ursa_stock.csv \
  --candidates data/candidates.json

pytest
ruff check src tests
```
