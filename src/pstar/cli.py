"""Command-line entry points."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from .chemistry import canonicalize
from .io import load_targets, read_jsonl, write_jsonl, write_portfolios
from .policies import FeatureTable, PaperFeatureProvider
from .portfolio import load_policy_specs, run_policy, run_portfolio
from .predictions import JsonPredictionSource
from .stock import Stock
from .validate import validate_route

DEFAULT_VALUE_MAP = Path(__file__).parent / "data" / "value_map.json"


def _resources(arguments):
    targets = [canonicalize(row) or row for row in load_targets(arguments.targets)]
    stock = Stock.from_csv(arguments.stock)
    predictions = JsonPredictionSource.from_file(arguments.candidates)
    features = (
        FeatureTable(arguments.descriptors)
        if arguments.descriptors
        else PaperFeatureProvider()
    )
    return targets, stock, predictions, features


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--targets", required=True, help="CSV with a smiles column")
    parser.add_argument("--stock", required=True, help="CSV with a smiles column")
    parser.add_argument("--candidates", required=True, help="product-to-candidates JSON")
    parser.add_argument("--descriptors", help="optional exact nine-feature JSON table")
    parser.add_argument("--value-map", default=str(DEFAULT_VALUE_MAP))
    parser.add_argument("--output", required=True)


def search_main() -> None:
    parser = argparse.ArgumentParser(description="Run one explicit P* policy")
    _common(parser)
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True, help="policy tag in the YAML config")
    args = parser.parse_args()
    document = yaml.safe_load(Path(args.config).read_text())
    specs = {row.tag: row for row in load_policy_specs(args.config)}
    spec = specs[args.policy]
    targets, stock, predictions, features = _resources(args)
    rows = []
    for index, target in enumerate(targets):
        result = run_policy(
            target,
            spec,
            predictions,
            stock,
            features=features,
            value_map=args.value_map,
            budget=int(document.get("budget", 500)),
            max_depth=int(document.get("max_depth", 30)),
        )
        row = result.to_dict()
        row.update({"idx": index, "target": target, "policy": spec.tag})
        rows.append(row)
    write_jsonl(args.output, rows)


def portfolio_main() -> None:
    parser = argparse.ArgumentParser(description="Run the frozen P* portfolio")
    _common(parser)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    document = yaml.safe_load(Path(args.config).read_text())
    targets, stock, predictions, features = _resources(args)
    portfolios = run_portfolio(
        targets,
        predictions,
        stock,
        features=features,
        value_map=args.value_map,
        policies=load_policy_specs(args.config),
        budget=int(document.get("budget", 500)),
        max_depth=int(document.get("max_depth", 30)),
        max_routes=int(document.get("max_routes", 10)),
    )
    write_portfolios(args.output, portfolios)


def validate_main() -> None:
    parser = argparse.ArgumentParser(description="Independently validate search routes")
    parser.add_argument("--results", required=True)
    parser.add_argument("--stock", required=True)
    parser.add_argument("--candidates", required=True)
    args = parser.parse_args()
    stock = Stock.from_csv(args.stock)
    predictions = JsonPredictionSource.from_file(args.candidates)
    failures = []
    checked = 0
    for row in read_jsonl(args.results):
        if "routes" in row:
            routes = [entry["route"] for entry in row["routes"]]
        elif row.get("route"):
            routes = [row["route"]]
        else:
            routes = []
        for route in routes:
            checked += 1
            result = validate_route(row["target"], route, predictions, stock)
            if not result.valid:
                failures.append({"idx": row.get("idx"), "reason": result.reason})
    print(json.dumps({"checked": checked, "invalid": len(failures), "failures": failures}, indent=2))
    if failures:
        raise SystemExit(1)


def analyze_main() -> None:
    parser = argparse.ArgumentParser(description="Summarize one or more result JSONL files")
    parser.add_argument("results", nargs="+")
    args = parser.parse_args()
    for path in args.results:
        rows = read_jsonl(path)
        solved = sum(bool(row.get("solved")) for row in rows)
        budgets = [row["solved_at"] for row in rows if row.get("solved_at") is not None]
        summary = {
            "path": path,
            "targets": len(rows),
            "solved": solved,
            "success": solved / len(rows) if rows else 0.0,
            "mean_solve_budget": sum(budgets) / len(budgets) if budgets else None,
        }
        print(json.dumps(summary, sort_keys=True))
