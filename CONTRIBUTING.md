# Contributing to pace-earthcare-matchups

Thank you for your interest in contributing to `pace-earthcare-matchups` (PEM)! Any help is greatly appreciated.

## Environment

First, you'll need some optional dependencies. Please use the following pip install command, rather than the one in the README:
```
pip install -e ".[maintainers]"
```

If you plan on modifying the documentation, you'll also need the `docs` dependencies:
```
pip install -e ".[maintainers,docs]"
```

## Formatting and Checking

If you've correctly installed the environment, then every time you call `git commit`, you should see that [`pre-commit`](https://pre-commit.com/) automatically runs the [`ruff`](https://docs.astral.sh/ruff/) formatter and checker. This should ensure you pass CI.

### Optional Type Checking

Static type checking with [`pyright`](https://github.com/microsoft/pyright) is configured in `pyproject.toml`. While this is not currently enforced in CI, type-checked code is preferred and greatly appreciated.

## Notebooks

Jupyter notebooks live in the [notebooks](notebooks) folder. Since `.ipynb` files don't play nicely with `git`, notebooks are stored as [MyST markdown](https://mystmd.org/). You can convert between notebooks and markdown as follows:
```
jupytext --to notebook notebooks/tutorial.md
jupytext --to myst notebooks/tutorial.ipynb
```

## AI Policy

The usage of AI tools to help write or review code is permissible for contributions to PEM. However, it is a requirement that all contributions are reviewed in their entirety by a human being. This means you should fully understand every line of code you contribute, and such code should be easily understandable by the maintainers of PEM. Fully autonomous contributions (e.g. a PR or entirely authored by an agent) are not allowed.

If you do use AI to help with your contribution to PEM, we ask that you disclose the nature of such usage in the pull request.
