# Contributing

Bug reports and pull requests are welcome.

Before opening a pull request:

```bash
ruff check src tests
pytest
```

Changes to the frozen policy list or benchmark protocol must be placed in a new
configuration or protocol document. Do not silently alter
`configs/pstar.yaml`.

Please do not commit model weights, commercial catalogues, reaction databases,
or benchmark data whose license does not allow redistribution.
