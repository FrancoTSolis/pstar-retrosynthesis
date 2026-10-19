"""Target and result file helpers."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

from .portfolio import TargetPortfolio


def load_targets(path: str | Path, smiles_column: str = "smiles") -> list[str]:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return []
    if smiles_column not in rows[0]:
        aliases = {key.lower(): key for key in rows[0]}
        smiles_column = aliases.get("smiles", aliases.get("target", smiles_column))
    return [row[smiles_column] for row in rows]


def write_jsonl(path: str | Path, rows: Iterable[dict]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def read_jsonl(path: str | Path) -> list[dict]:
    with open(path) as handle:
        return [json.loads(line) for line in handle if line.strip()]


def write_portfolios(output: str | Path, portfolios: list[TargetPortfolio]) -> None:
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    for policy in {run.policy for target in portfolios for run in target.policies}:
        rows = []
        for index, target in enumerate(portfolios):
            run = next(row for row in target.policies if row.policy == policy)
            row = run.result.to_dict()
            row.update({"idx": index, "target": target.target, "policy": policy})
            rows.append(row)
        write_jsonl(root / f"{policy}.jsonl", rows)

    portfolio_rows = []
    for index, target in enumerate(portfolios):
        portfolio_rows.append(
            {
                "idx": index,
                "target": target.target,
                "n_unique_routes": len(target.unique_routes),
                "routes": [
                    {"policy": policy, "route": route}
                    for policy, route in target.unique_routes
                ],
            }
        )
    write_jsonl(root / "portfolio.jsonl", portfolio_rows)
