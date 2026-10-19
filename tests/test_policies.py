from pstar.candidate import Candidate
from pstar.engine import SearchContext
from pstar.policies import building_block_scorer, prior_gated_building_block_scorer
from pstar.stock import Stock


def test_building_block_rule_prefers_more_stock_precursors():
    stock = Stock.from_smiles(["C", "CO"])
    context = SearchContext(stock=stock, target="CCO")
    candidates = (
        Candidate(("C", "CO"), 0.01),
        Candidate(("CC",), 0.99),
    )
    scores = building_block_scorer()(context.target, candidates, context)
    assert scores[0] > scores[1]


def test_prior_gate_rejects_low_probability_bb_candidate():
    stock = Stock.from_smiles(["C", "CO"])
    context = SearchContext(stock=stock, target="CCO")
    candidates = (
        Candidate(("C", "CO"), 0.001),
        Candidate(("CC",), 0.5),
    )
    scores = prior_gated_building_block_scorer(threshold=0.02)(
        context.target,
        candidates,
        context,
    )
    assert scores[1] > scores[0]
