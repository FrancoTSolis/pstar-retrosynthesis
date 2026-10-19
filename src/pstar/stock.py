"""Canonical-SMILES building-block catalogue."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .chemistry import canonicalize


@dataclass(slots=True)
class Stock:
    smiles: frozenset[str]

    def __contains__(self, smiles: str) -> bool:
        canon = canonicalize(smiles)
        return canon in self.smiles if canon else False

    def __call__(self, smiles: str) -> bool:
        return smiles in self

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        smiles_column: str = "smiles",
        canonicalize_rows: bool = True,
    ) -> "Stock":
        values: set[str] = set()
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                raw = row[smiles_column]
                value = canonicalize(raw) if canonicalize_rows else raw
                if value:
                    values.add(value)
        return cls(frozenset(values))

    @classmethod
    def from_smiles(cls, rows: list[str] | tuple[str, ...]) -> "Stock":
        values = {canonicalize(row) for row in rows}
        return cls(frozenset(value for value in values if value))
