"""One-step retrosynthesis candidate representation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class Candidate:
    """A product decomposition proposed by a single-step model."""

    reactants: tuple[str, ...]
    score: float
    template: str | None = None

    @classmethod
    def from_dict(cls, row: dict[str, Any]) -> "Candidate":
        return cls(
            reactants=tuple(sorted(str(x) for x in row["reactants"])),
            score=float(row.get("score", 0.0)),
            template=row.get("template"),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "reactants": list(self.reactants),
            "score": self.score,
            "template": self.template,
        }
