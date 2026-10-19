"""Independent validation of emitted routes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .portfolio import route_signature
from .predictions import PredictionSource


@dataclass(frozen=True, slots=True)
class ValidationResult:
    valid: bool
    reason: str | None
    depth: int | None


def validate_route(
    target: str,
    route: Sequence[tuple[str, Sequence[str]]],
    predictions: PredictionSource,
    stock,
    *,
    max_depth: int = 30,
) -> ValidationResult:
    steps: dict[str, tuple[str, ...]] = {}
    for product, reactants in route:
        if product in steps:
            return ValidationResult(False, f"duplicate product step: {product}", None)
        steps[product] = tuple(sorted(reactants))
        candidate_sets = {
            tuple(sorted(candidate.reactants))
            for candidate in predictions.get(product)
        }
        if steps[product] not in candidate_sets:
            return ValidationResult(
                False,
                f"step is not a candidate for product: {product}",
                None,
            )

    if target not in steps and not stock(target):
        return ValidationResult(False, "target has no route step", None)

    leaves = {
        reactant
        for reactants in steps.values()
        for reactant in reactants
        if reactant not in steps
    }
    missing = sorted(row for row in leaves if not stock(row))
    if missing:
        return ValidationResult(False, f"{len(missing)} route leaves are not in stock", None)

    memo: dict[str, int] = {}

    def depth(smiles: str, ancestors: frozenset[str]) -> int:
        if smiles not in steps:
            return 0
        if smiles in ancestors:
            raise ValueError(f"cycle at {smiles}")
        if smiles not in memo:
            memo[smiles] = 1 + max(
                depth(row, ancestors | {smiles})
                for row in steps[smiles]
            )
        return memo[smiles]

    try:
        route_depth = depth(target, frozenset())
    except ValueError as error:
        return ValidationResult(False, str(error), None)
    if route_depth > max_depth:
        return ValidationResult(False, f"route depth {route_depth} exceeds {max_depth}", route_depth)
    route_signature(route)
    return ValidationResult(True, None, route_depth)
