"""AND-OR best-first retrosynthesis search with pluggable local policies."""

from __future__ import annotations

import math
from collections import deque
from dataclasses import asdict, dataclass
from typing import Callable, Sequence

from .candidate import Candidate
from .predictions import PredictionSource

INF = float("inf")
ReactionScorer = Callable[[str, Sequence[Candidate], "SearchContext"], list[float]]
LeafValue = Callable[[Sequence[str], "SearchContext"], list[float]]
Width = Callable[[Sequence[float]], int]
StockCheck = Callable[[str], bool]


@dataclass(slots=True)
class SearchContext:
    stock: StockCheck
    target: str


@dataclass(slots=True)
class SearchResult:
    solved: bool
    solved_at: int | None
    budget: int
    route_len: int | None
    route_depth: int | None
    route: list[tuple[str, list[str]]] | None
    stats: dict[str, int]
    n_molecules: int

    def to_dict(self) -> dict:
        return asdict(self)


def _log_softmax(values: Sequence[float]) -> list[float]:
    maximum = max(values)
    log_z = maximum + math.log(sum(math.exp(value - maximum) for value in values))
    return [value - log_z for value in values]


class _Molecule:
    def __init__(self, smiles: str, index: int, building_block: bool, depth: int) -> None:
        self.smiles = smiles
        self.index = index
        self.building_block = building_block
        self.depth = depth
        self.parents: list[_Reaction] = []
        self.children: list[_Reaction] = []
        self.expanded = False
        self.dead = False
        self.value = 0.0
        self.route_number = 0.0
        self.reserve: list[tuple[float, Candidate]] = []
        self.solved = building_block


class _Reaction:
    def __init__(
        self,
        parent: _Molecule,
        children: list[_Molecule],
        cost: float,
        candidate: Candidate,
    ) -> None:
        self.parent = parent
        self.children = children
        self.cost = cost
        self.candidate = candidate
        self.route_number = INF
        self.solved = False

    def update(self) -> None:
        self.route_number = self.cost + sum(child.route_number for child in self.children)
        self.solved = all(child.solved for child in self.children)


class Search:
    """Run one deterministic search policy for one target."""

    def __init__(
        self,
        target: str,
        predictions: PredictionSource,
        stock: StockCheck,
        scorer: ReactionScorer,
        value: LeafValue,
        width: Width,
        *,
        max_budget: int = 500,
        max_depth: int = 30,
    ) -> None:
        self.predictions = predictions
        self.stock = stock
        self.scorer = scorer
        self.value = value
        self.width = width
        self.max_budget = max_budget
        self.max_depth = max_depth
        self.context = SearchContext(stock=stock, target=target)
        self.molecules: dict[str, _Molecule] = {}
        self.budget = 0
        self.solved_at: int | None = None
        self.stats = {
            "cycle_rejected": 0,
            "dead_leaves": 0,
            "reserve_opened": 0,
            "candidates_seen": 0,
            "reactions_opened": 0,
            "expansions": 0,
        }
        self.root = self._molecule(target, 0)
        self._initialize_values([self.root])
        if self.root.solved:
            self.solved_at = 0

    def _molecule(self, smiles: str, depth: int) -> _Molecule:
        existing = self.molecules.get(smiles)
        if existing is not None:
            if depth < existing.depth:
                self._lower_depth(existing, depth)
            return existing
        node = _Molecule(smiles, len(self.molecules), bool(self.stock(smiles)), depth)
        self.molecules[smiles] = node
        return node

    def _lower_depth(self, start: _Molecule, depth: int) -> None:
        queue = deque([(start, depth)])
        while queue:
            node, new_depth = queue.popleft()
            if new_depth >= node.depth:
                continue
            node.depth = new_depth
            for reaction in node.children:
                for child in reaction.children:
                    if new_depth + 1 < child.depth:
                        queue.append((child, new_depth + 1))

    def _is_ancestor_or_self(self, ancestor: _Molecule, node: _Molecule) -> bool:
        if ancestor is node:
            return True
        seen = {node.index}
        stack = [node]
        while stack:
            current = stack.pop()
            for parent_reaction in current.parents:
                parent = parent_reaction.parent
                if parent is ancestor:
                    return True
                if parent.index not in seen:
                    seen.add(parent.index)
                    stack.append(parent)
        return False

    def _creates_cycle(self, product: _Molecule, candidate: Candidate) -> bool:
        return any(
            (existing := self.molecules.get(smiles)) is not None
            and self._is_ancestor_or_self(existing, product)
            for smiles in candidate.reactants
        )

    def _initialize_values(self, molecules: Sequence[_Molecule]) -> None:
        open_rows = [row for row in molecules if not row.building_block and not row.expanded]
        if open_rows:
            values = self.value([row.smiles for row in open_rows], self.context)
            if len(values) != len(open_rows):
                raise ValueError("leaf-value function returned the wrong number of values")
            for row, value in zip(open_rows, values):
                row.value = max(0.0, float(value))
        for row in molecules:
            row.route_number = self._leaf_route_number(row)

    def _leaf_route_number(self, molecule: _Molecule) -> float:
        if molecule.building_block:
            return 0.0
        if molecule.dead:
            return INF
        if not molecule.expanded:
            return INF if molecule.depth >= self.max_depth else molecule.value
        return min((reaction.route_number for reaction in molecule.children), default=INF)

    def _open(self, product: _Molecule, items: Sequence[tuple[float, Candidate]]) -> None:
        new_molecules: list[_Molecule] = []
        reactions: list[_Reaction] = []
        for cost, candidate in items:
            self.stats["candidates_seen"] += 1
            if self._creates_cycle(product, candidate):
                self.stats["cycle_rejected"] += 1
                continue
            old_count = len(self.molecules)
            children = [self._molecule(smiles, product.depth + 1) for smiles in candidate.reactants]
            for child in children:
                if child.index >= old_count and child not in new_molecules:
                    new_molecules.append(child)
            reaction = _Reaction(product, children, cost, candidate)
            for child in children:
                child.parents.append(reaction)
            product.children.append(reaction)
            reactions.append(reaction)
            self.stats["reactions_opened"] += 1
        self._initialize_values(new_molecules)
        for reaction in reactions:
            reaction.update()

    def _refill(self, molecule: _Molecule) -> None:
        while molecule.reserve and all(row.route_number == INF for row in molecule.children):
            item = molecule.reserve.pop(0)
            self.stats["reserve_opened"] += 1
            self._open(molecule, [item])

    def _expand(self, molecule: _Molecule) -> None:
        candidates = list(self.predictions.get(molecule.smiles))
        self.budget += 1
        self.stats["expansions"] += 1
        molecule.expanded = True

        valid: list[Candidate] = []
        for candidate in candidates:
            if self._creates_cycle(molecule, candidate):
                self.stats["cycle_rejected"] += 1
            else:
                valid.append(candidate)
        keep = getattr(self.predictions, "keep", None)
        if keep is not None:
            valid = valid[: int(keep)]
        if not valid:
            molecule.dead = True
            self.stats["dead_leaves"] += 1
            return

        logits = self.scorer(molecule.smiles, valid, self.context)
        if len(logits) != len(valid):
            raise ValueError("reaction scorer returned the wrong number of values")
        log_probabilities = _log_softmax(logits)
        order = sorted(range(len(valid)), key=lambda index: -logits[index])
        items = [(-log_probabilities[index], valid[index]) for index in order]
        count = max(1, min(len(items), int(self.width([logits[index] for index in order]))))
        molecule.reserve = items[count:]
        self._open(molecule, items[:count])
        self._refill(molecule)
        if not molecule.children:
            molecule.dead = True
            self.stats["dead_leaves"] += 1

    def _recompute(self, start: _Molecule) -> None:
        queue = deque([start])
        queued = {start.index}
        first = True
        while queue:
            molecule = queue.popleft()
            queued.discard(molecule.index)
            if molecule.expanded and not molecule.building_block:
                for reaction in molecule.children:
                    reaction.update()
                self._refill(molecule)
                route_number = min(
                    (reaction.route_number for reaction in molecule.children),
                    default=INF,
                )
                solved = any(reaction.solved for reaction in molecule.children)
            else:
                route_number = self._leaf_route_number(molecule)
                solved = molecule.building_block
            if first or route_number != molecule.route_number or solved != molecule.solved:
                first = False
                molecule.route_number = route_number
                molecule.solved = solved
                for reaction in molecule.parents:
                    parent = reaction.parent
                    if parent.index not in queued:
                        queued.add(parent.index)
                        queue.append(parent)

    def _select(self) -> _Molecule | None:
        best: _Molecule | None = None
        seen: set[int] = set()
        stack = [self.root]
        while stack:
            molecule = stack.pop()
            if molecule.index in seen or molecule.solved:
                continue
            seen.add(molecule.index)
            if not molecule.expanded:
                if molecule.route_number < INF and (
                    best is None or molecule.index < best.index
                ):
                    best = molecule
                continue
            live = [row for row in molecule.children if row.route_number < INF]
            if live:
                stack.extend(min(live, key=lambda row: row.route_number).children)
        return best

    def run(self) -> SearchResult:
        while not self.root.solved:
            if self.budget >= self.max_budget or self.root.route_number == INF:
                break
            molecule = self._select()
            if molecule is None:
                break
            self._expand(molecule)
            self._recompute(molecule)
            if self.root.solved:
                self.solved_at = self.budget
        return self.result()

    def route(self) -> tuple[list[tuple[str, list[str]]], int] | None:
        if not self.root.solved:
            return None
        memo: dict[str, float] = {}

        def cost(molecule: _Molecule) -> float:
            if molecule.building_block:
                return 0.0
            if molecule.smiles not in memo:
                memo[molecule.smiles] = min(
                    (
                        reaction.cost + sum(cost(child) for child in reaction.children)
                        for reaction in molecule.children
                        if reaction.solved
                    ),
                    default=INF,
                )
            return memo[molecule.smiles]

        steps: list[tuple[str, list[str]]] = []
        seen: set[str] = set()
        depths: dict[str, int] = {}

        def visit(molecule: _Molecule) -> int:
            if molecule.building_block:
                return 0
            if molecule.smiles in depths:
                return depths[molecule.smiles]
            reaction = min(
                (row for row in molecule.children if row.solved),
                key=lambda row: row.cost + sum(cost(child) for child in row.children),
            )
            if molecule.smiles not in seen:
                seen.add(molecule.smiles)
                steps.append((molecule.smiles, [child.smiles for child in reaction.children]))
            depths[molecule.smiles] = 1 + max(visit(child) for child in reaction.children)
            return depths[molecule.smiles]

        return steps, visit(self.root)

    def result(self) -> SearchResult:
        route_result = self.route()
        route, depth = route_result if route_result else (None, None)
        return SearchResult(
            solved=self.root.solved,
            solved_at=self.solved_at,
            budget=self.budget,
            route_len=len(route) if route else None,
            route_depth=depth,
            route=route,
            stats=dict(self.stats),
            n_molecules=len(self.molecules),
        )
