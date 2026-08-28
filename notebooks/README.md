# Notebooks

Intentionally empty. Every analysis in the paper is a command-line entry point under `src/` or
`scripts/`, so that results are reproducible without notebook state.

Start from `scripts/reproduce_all.sh --from-precomputed`, or call a module directly:

```bash
cd src && python3 -m analysis.evaluator_agreement --help
```
