from pstar.candidate import Candidate
from pstar.engine import Search
from pstar.policies import prior_scorer, width_all, zero_value
from pstar.predictions import JsonPredictionSource
from pstar.stock import Stock


def test_search_stops_at_first_complete_route():
    source = JsonPredictionSource(
        {
            "CCO": [
                {"reactants": ["CC"], "score": 0.8},
                {"reactants": ["CO"], "score": 0.2},
            ],
            "CC": [{"reactants": ["C"], "score": 1.0}],
        }
    )
    stock = Stock.from_smiles(["C", "CO"])
    result = Search(
        "CCO",
        source,
        stock,
        prior_scorer,
        zero_value,
        width_all,
        max_budget=10,
    ).run()
    assert result.solved
    assert result.solved_at == 1
    assert result.route == [("CCO", ["CO"])]


def test_cycle_is_rejected():
    source = JsonPredictionSource(
        {
            "CCO": [{"reactants": ["CC"], "score": 1.0}],
            "CC": [
                {"reactants": ["CCO"], "score": 0.9},
                {"reactants": ["C"], "score": 0.1},
            ],
        }
    )
    stock = Stock.from_smiles(["C"])
    result = Search(
        "CCO",
        source,
        stock,
        prior_scorer,
        zero_value,
        width_all,
        max_budget=10,
    ).run()
    assert result.solved
    assert result.stats["cycle_rejected"] == 1


def test_budget_counts_distinct_revealed_molecules():
    source = JsonPredictionSource({"CCO": [{"reactants": ["CC"], "score": 1.0}]})
    result = Search(
        "CCO",
        source,
        Stock.from_smiles([]),
        prior_scorer,
        zero_value,
        width_all,
        max_budget=1,
    ).run()
    assert not result.solved
    assert result.budget == 1


def test_prior_scorer_returns_log_probability():
    values = prior_scorer(
        "CCO",
        (Candidate(("CC",), 0.75), Candidate(("CO",), 0.25)),
        None,
    )
    assert values[0] > values[1]
