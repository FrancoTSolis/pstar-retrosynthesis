"""Reaction scorers, leaf values, and frozen policy building blocks."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Callable, Protocol, Sequence

from rdkit import Chem, RDConfig

from .candidate import Candidate
from .chemistry import basic_descriptors, heavy_atoms
from .engine import LeafValue, ReactionScorer, SearchContext, Width


def prior_scorer(
    product: str,
    candidates: Sequence[Candidate],
    context: SearchContext,
) -> list[float]:
    del product, context
    return [math.log(max(candidate.score, 1e-12)) for candidate in candidates]


def building_block_scorer(
    purchasable_weight: float = 2.0,
    synthesis_weight: float = 1.0,
    retention_weight: float = 5.0,
) -> ReactionScorer:
    def score(
        product: str,
        candidates: Sequence[Candidate],
        context: SearchContext,
    ) -> list[float]:
        product_size = max(1, heavy_atoms(product))
        values: list[float] = []
        for candidate in candidates:
            open_rows = [row for row in candidate.reactants if not context.stock(row)]
            purchasable = len(candidate.reactants) - len(open_rows)
            retention = max((heavy_atoms(row) for row in open_rows), default=0) / product_size
            values.append(
                purchasable_weight * purchasable
                - synthesis_weight * len(open_rows)
                - retention_weight * retention
            )
        return values

    return score


def prior_gated_building_block_scorer(
    *,
    threshold: float = 0.02,
    purchasable_weight: float = 2.0,
    synthesis_weight: float = 1.0,
    retention_weight: float = 5.0,
) -> ReactionScorer:
    inner = building_block_scorer(
        purchasable_weight,
        synthesis_weight,
        retention_weight,
    )

    def score(
        product: str,
        candidates: Sequence[Candidate],
        context: SearchContext,
    ) -> list[float]:
        values = inner(product, candidates, context)
        probabilities = [candidate.score for candidate in candidates]
        if not any(probability >= threshold for probability in probabilities):
            return [math.log(max(probability, 1e-12)) for probability in probabilities]
        floor = min(values) - 1.0
        return [
            value
            if probability >= threshold
            else floor - 1.0 + 1e-3 * math.log(max(probability, 1e-12))
            for value, probability in zip(values, probabilities)
        ]

    return score


def zero_value(smiles: Sequence[str], context: SearchContext) -> list[float]:
    del context
    return [0.0] * len(smiles)


def constant_value(value: float) -> LeafValue:
    def evaluate(smiles: Sequence[str], context: SearchContext) -> list[float]:
        del context
        return [float(value)] * len(smiles)

    return evaluate


def external_value(predict: Callable[[Sequence[str]], Sequence[float]], scale: float) -> LeafValue:
    def evaluate(smiles: Sequence[str], context: SearchContext) -> list[float]:
        del context
        return [scale * max(0.0, float(row)) for row in predict(smiles)]

    return evaluate


class FeatureProvider(Protocol):
    def __call__(self, smiles: str) -> Sequence[float]: ...


class FeatureTable:
    """Nine-descriptor feature table keyed by canonical SMILES."""

    def __init__(self, path: str | Path) -> None:
        with open(path) as handle:
            self.rows = json.load(handle)

    def __call__(self, smiles: str) -> Sequence[float]:
        try:
            return self.rows[smiles]
        except KeyError as error:
            raise KeyError(f"no descriptor row for {smiles}") from error


class PaperFeatureProvider:
    """Compute the paper's nine features.

    ``scscore`` and ``desp_value`` callbacks are required for exact benchmark
    reproduction. Defaults keep the toy and unit-test workflows self-contained.
    """

    def __init__(
        self,
        *,
        scscore: Callable[[str], float] | None = None,
        desp_value: Callable[[str], float] | None = None,
    ) -> None:
        self.scscore = scscore or (lambda _: 3.0)
        self.desp_value = desp_value or (lambda _: 0.0)
        sys.path.append(str(Path(RDConfig.RDContribDir) / "SA_Score"))
        try:
            import sascorer  # type: ignore
        except ImportError:
            from rdkit.Contrib.SA_Score import sascorer  # type: ignore
        self.sascorer = sascorer

    def __call__(self, smiles: str) -> Sequence[float]:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return (5.0, 3.0, 2.0, 2.0, 2.0, 2.0, 5.0, 3.0, 3.0)
        base = list(basic_descriptors(smiles))
        return (
            *base,
            float(self.scscore(smiles)),
            float(self.sascorer.calculateScore(molecule)) / 3.0,
            max(0.0, float(self.desp_value(smiles))) / 3.0,
        )


class LogisticValue:
    def __init__(
        self,
        map_path: str | Path,
        feature_provider: FeatureProvider,
        alpha: float,
    ) -> None:
        with open(map_path) as handle:
            model = json.load(handle)
        self.weights = tuple(float(row) for row in model["w"])
        self.bias = float(model["b"])
        self.mean = tuple(float(row) for row in model["mu"])
        self.std = tuple(float(row) for row in model["sd"])
        self.features = feature_provider
        self.alpha = float(alpha)

    def __call__(self, smiles: Sequence[str], context: SearchContext) -> list[float]:
        del context
        output: list[float] = []
        for row in smiles:
            features = self.features(row)
            standardized = [
                (value - mean) / (std or 1.0)
                for value, mean, std in zip(features, self.mean, self.std)
            ]
            logit = self.bias + sum(
                weight * value for weight, value in zip(self.weights, standardized)
            )
            probability = 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, logit))))
            output.append(self.alpha * -math.log(max(1e-6, probability)))
        return output


def width_all(logits: Sequence[float]) -> int:
    return len(logits)


def width_top_k(count: int) -> Width:
    return lambda logits: min(count, len(logits))
