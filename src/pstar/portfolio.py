"""Frozen policy definitions and portfolio execution."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import yaml

from .engine import LeafValue, ReactionScorer, Search, SearchResult
from .policies import (
    FeatureProvider,
    LogisticValue,
    building_block_scorer,
    constant_value,
    prior_gated_building_block_scorer,
    prior_scorer,
    width_all,
)
from .predictions import MatchedView, PredictionSource


@dataclass(frozen=True, slots=True)
class PolicySpec:
    tag: str
    reaction_score: str
    leaf_value: str
    alpha: float | None = None
    view: str = "full"
    seed: int = 0
    threshold: float = 0.02

    @classmethod
    def from_dict(cls, row: dict) -> "PolicySpec":
        return cls(
            tag=str(row["tag"]),
            reaction_score=str(row["reaction_score"]),
            leaf_value=str(row["leaf_value"]),
            alpha=float(row["alpha"]) if row.get("alpha") is not None else None,
            view=str(row.get("view", "full")),
            seed=int(row.get("seed", 0)),
            threshold=float(row.get("threshold", 0.02)),
        )


def frozen_policies() -> tuple[PolicySpec, ...]:
    return (
        PolicySpec("P01_prior_L16", "prior", "learned", 16),
        PolicySpec("P02_prior_const2", "prior", "constant2"),
        PolicySpec("P03_bb0_L64", "bb_no_reward", "learned", 64),
        PolicySpec("P04_bb0_L8", "bb_no_reward", "learned", 8),
        PolicySpec("P05_prior_d4", "prior", "desp"),
        PolicySpec("P06_bb_L8", "bb", "learned", 8),
        PolicySpec("P07_bbg002_L16", "prior_gated_bb", "learned", 16, threshold=0.02),
        PolicySpec("P08_Cstar", "bb", "learned", 16),
        PolicySpec("P09_Mstar", "bb", "learned", 8, view="matched", seed=0),
        PolicySpec("P10_prior_scs1", "prior", "scscore"),
    )


def load_policy_specs(path: str | Path) -> tuple[PolicySpec, ...]:
    with open(path) as handle:
        document = yaml.safe_load(handle)
    return tuple(PolicySpec.from_dict(row) for row in document["policies"])


def _reaction_scorer(spec: PolicySpec) -> ReactionScorer:
    if spec.reaction_score == "prior":
        return prior_scorer
    if spec.reaction_score == "bb":
        return building_block_scorer(2.0, 1.0, 5.0)
    if spec.reaction_score == "bb_no_reward":
        return building_block_scorer(0.0, 1.0, 5.0)
    if spec.reaction_score == "prior_gated_bb":
        return prior_gated_building_block_scorer(
            threshold=spec.threshold,
            purchasable_weight=2.0,
            synthesis_weight=1.0,
            retention_weight=5.0,
        )
    raise ValueError(f"unknown reaction score: {spec.reaction_score}")


def _leaf_value(
    spec: PolicySpec,
    value_map: str | Path,
    features: FeatureProvider,
) -> LeafValue:
    if spec.leaf_value == "learned":
        if spec.alpha is None:
            raise ValueError(f"{spec.tag} requires alpha")
        return LogisticValue(value_map, features, spec.alpha)
    if spec.leaf_value == "constant2":
        return constant_value(2.0)
    if spec.leaf_value == "desp":
        def value(smiles, context):
            del context
            return [4.0 * max(0.0, float(features(row)[8]) * 3.0) for row in smiles]
        return value
    if spec.leaf_value == "scscore":
        def value(smiles, context):
            del context
            return [max(0.0, float(features(row)[6]) - 1.0) for row in smiles]
        return value
    raise ValueError(f"unknown leaf value: {spec.leaf_value}")


def route_signature(route: Sequence[tuple[str, Sequence[str]]]) -> tuple:
    return tuple(
        sorted((product, tuple(sorted(reactants))) for product, reactants in route)
    )


@dataclass(slots=True)
class PolicyRun:
    policy: str
    result: SearchResult


@dataclass(slots=True)
class TargetPortfolio:
    target: str
    policies: list[PolicyRun]
    unique_routes: list[tuple[str, list[tuple[str, list[str]]]]]


def run_policy(
    target: str,
    spec: PolicySpec,
    predictions: PredictionSource,
    stock: Callable[[str], bool],
    *,
    features: FeatureProvider,
    value_map: str | Path,
    budget: int = 500,
    max_depth: int = 30,
) -> SearchResult:
    source = (
        MatchedView(predictions, seed=spec.seed, keep=20)
        if spec.view == "matched"
        else predictions
    )
    search = Search(
        target,
        source,
        stock,
        _reaction_scorer(spec),
        _leaf_value(spec, value_map, features),
        width_all,
        max_budget=budget,
        max_depth=max_depth,
    )
    return search.run()


def run_portfolio(
    targets: Iterable[str],
    predictions: PredictionSource,
    stock: Callable[[str], bool],
    *,
    features: FeatureProvider,
    value_map: str | Path,
    policies: Sequence[PolicySpec] | None = None,
    budget: int = 500,
    max_depth: int = 30,
    max_routes: int = 10,
) -> list[TargetPortfolio]:
    specs = tuple(policies or frozen_policies())
    portfolios: list[TargetPortfolio] = []
    for target in targets:
        runs: list[PolicyRun] = []
        unique: list[tuple[str, list[tuple[str, list[str]]]]] = []
        signatures: set[tuple] = set()
        for spec in specs:
            result = run_policy(
                target,
                spec,
                predictions,
                stock,
                features=features,
                value_map=value_map,
                budget=budget,
                max_depth=max_depth,
            )
            runs.append(PolicyRun(spec.tag, result))
            if result.route and len(unique) < max_routes:
                signature = route_signature(result.route)
                if signature not in signatures:
                    signatures.add(signature)
                    unique.append((spec.tag, result.route))
        portfolios.append(TargetPortfolio(target, runs, unique))
    return portfolios
