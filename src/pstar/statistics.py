"""Small dependency-free statistical utilities used in the paper."""

from __future__ import annotations

import math
import random
from collections.abc import Sequence


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        raise ValueError("total must be positive")
    probability = successes / total
    denominator = 1.0 + z * z / total
    center = probability + z * z / (2.0 * total)
    radius = z * math.sqrt(
        probability * (1.0 - probability) / total + z * z / (4.0 * total * total)
    )
    return (center - radius) / denominator, (center + radius) / denominator


def exact_mcnemar(wins: int, losses: int) -> float:
    """Exact two-sided binomial McNemar p-value."""
    discordant = wins + losses
    if discordant == 0:
        return 1.0
    smaller = min(wins, losses)
    probability = sum(
        math.comb(discordant, value) * 0.5**discordant
        for value in range(smaller + 1)
    )
    return min(1.0, 2.0 * probability)


def sign_permutation_p(
    values: Sequence[float],
    *,
    samples: int = 200_000,
    seed: int = 0,
) -> float:
    rng = random.Random(seed)
    observed = abs(sum(values))
    exceed = 0
    for _ in range(samples):
        statistic = sum(value if rng.random() < 0.5 else -value for value in values)
        exceed += abs(statistic) >= observed - 1e-12
    return (exceed + 1) / (samples + 1)


def bootstrap_mean_interval(
    values: Sequence[float],
    *,
    samples: int = 20_000,
    seed: int = 0,
) -> tuple[float, float]:
    if not values:
        raise ValueError("values cannot be empty")
    rng = random.Random(seed)
    count = len(values)
    estimates = sorted(
        sum(values[rng.randrange(count)] for _ in range(count)) / count
        for _ in range(samples)
    )
    return estimates[int(0.025 * samples)], estimates[int(0.975 * samples)]


def holm_adjust(p_values: Sequence[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=lambda index: p_values[index])
    output = [1.0] * len(p_values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, p_values[index] * (len(order) - rank)))
        output[index] = running
    return output
