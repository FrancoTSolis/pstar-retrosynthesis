"""Export policy JSONL files to a RetroCast collected-routes archive."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import read_jsonl
from .portfolio import route_signature


def main() -> None:
    parser = argparse.ArgumentParser(description="Export P* routes for official URSA")
    parser.add_argument("--results", required=True, nargs="+", help="policy JSONL files in order")
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-routes", type=int, default=10)
    args = parser.parse_args()
    try:
        from retrocast import chem
        from retrocast.io import save_collected_routes
        from retrocast.models.route import Molecule, Reaction, Route
        from retrocast.typing import InChIKeyStr, SmilesStr
    except ImportError as error:
        raise ImportError('install the URSA extra with pip install -e ".[ursa]"') from error

    result_sets = [sorted(read_jsonl(path), key=lambda row: row["idx"]) for path in args.results]
    expected = [(row["idx"], row["target"]) for row in result_sets[0]]
    for path, rows in zip(args.results[1:], result_sets[1:]):
        if [(row["idx"], row["target"]) for row in rows] != expected:
            raise ValueError(f"target rows differ in {path}")

    by_target: dict[str, list[list]] = {}
    seen: dict[str, set[tuple]] = {}
    for rows in result_sets:
        for row in rows:
            route = row.get("route")
            if not row.get("solved") or not route:
                continue
            current = by_target.setdefault(row["target"], [])
            if len(current) >= args.max_routes:
                continue
            signature = route_signature(route)
            target_seen = seen.setdefault(row["target"], set())
            if signature not in target_seen:
                target_seen.add(signature)
                current.append(route)

    def molecule(smiles, steps, ancestors=frozenset()):
        canonical = chem.canonicalize_smiles(smiles)
        if canonical in ancestors:
            raise ValueError(f"cycle at {canonical}")
        reactants = steps.get(canonical)
        product_of = None
        if reactants is not None:
            product_of = Reaction(
                reactants=[
                    molecule(row, steps, ancestors | {canonical})
                    for row in reactants
                ]
            )
        return Molecule(
            smiles=SmilesStr(canonical),
            inchikey=InChIKeyStr(chem.get_inchi_key(canonical)),
            product_of=product_of,
        )

    collected = {}
    for target, routes in by_target.items():
        target_routes = []
        for route in routes:
            steps = {
                chem.canonicalize_smiles(product): [
                    chem.canonicalize_smiles(row) for row in reactants
                ]
                for product, reactants in route
            }
            target_routes.append(Route(target=molecule(target, steps)))
        if target_routes:
            collected[target] = target_routes
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    save_collected_routes(collected, args.output)
