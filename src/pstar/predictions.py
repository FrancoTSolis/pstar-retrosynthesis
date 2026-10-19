"""Prediction sources and deterministic matched-view sampling."""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path
from typing import Protocol, Sequence

from .candidate import Candidate
from .chemistry import canonicalize


class PredictionSource(Protocol):
    def get(self, smiles: str) -> Sequence[Candidate]: ...


class JsonPredictionSource:
    """In-memory product-to-candidate mapping loaded from JSON."""

    def __init__(self, rows: dict[str, list[dict]]) -> None:
        table: dict[str, tuple[Candidate, ...]] = {}
        for product, candidates in rows.items():
            canon = canonicalize(product) or product
            merged: dict[tuple[str, ...], Candidate] = {}
            for raw in candidates:
                candidate = Candidate.from_dict(raw)
                key = candidate.reactants
                if key in merged:
                    old = merged[key]
                    merged[key] = Candidate(key, old.score + candidate.score, old.template)
                else:
                    merged[key] = candidate
            total = sum(max(0.0, candidate.score) for candidate in merged.values())
            normalized = [
                Candidate(
                    candidate.reactants,
                    candidate.score / total if total > 0 else 1.0 / max(1, len(merged)),
                    candidate.template,
                )
                for candidate in merged.values()
            ]
            table[canon] = tuple(sorted(normalized, key=lambda row: -row.score))
        self.table = table

    @classmethod
    def from_file(cls, path: str | Path) -> "JsonPredictionSource":
        with open(path) as handle:
            return cls(json.load(handle))

    def get(self, smiles: str) -> tuple[Candidate, ...]:
        canon = canonicalize(smiles) or smiles
        return self.table.get(canon, ())


class MatchedView:
    """Seeded random ordering that hides scores and retains up to ``keep`` rows."""

    def __init__(self, source: PredictionSource, seed: int = 0, keep: int = 20) -> None:
        self.source = source
        self.seed = seed
        self.keep = keep

    def get(self, smiles: str) -> tuple[Candidate, ...]:
        rows = list(self.source.get(smiles))
        digest = hashlib.sha256(f"{self.seed}|{smiles}".encode()).hexdigest()
        random.Random(int(digest[:16], 16)).shuffle(rows)
        probability = 1.0 / max(1, len(rows))
        return tuple(Candidate(row.reactants, probability, row.template) for row in rows)
