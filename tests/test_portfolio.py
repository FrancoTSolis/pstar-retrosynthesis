from pathlib import Path

from pstar.io import load_targets
from pstar.policies import PaperFeatureProvider
from pstar.portfolio import run_portfolio
from pstar.predictions import JsonPredictionSource
from pstar.stock import Stock
from pstar.validate import validate_route


ROOT = Path(__file__).parents[1]


def test_frozen_portfolio_runs_on_toy_data():
    targets = load_targets(ROOT / "examples/toy/targets.csv")
    stock = Stock.from_csv(ROOT / "examples/toy/stock.csv")
    predictions = JsonPredictionSource.from_file(ROOT / "examples/toy/candidates.json")
    portfolios = run_portfolio(
        targets,
        predictions,
        stock,
        features=PaperFeatureProvider(),
        value_map=ROOT / "models/value_map.json",
        budget=20,
        max_depth=10,
    )
    assert len(portfolios) == 2
    assert all(len(row.policies) == 10 for row in portfolios)
    assert all(row.unique_routes for row in portfolios)
    for row in portfolios:
        for _, route in row.unique_routes:
            assert validate_route(row.target, route, predictions, stock).valid
