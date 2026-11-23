# Quantum-to-Einstein gravity research

An exploratory research notebook, not a demonstrated unified theory.

## Start here

1. [docs/README.md](docs/README.md) — documentation map  
2. [docs/OVERNIGHT_SUMMARY.md](docs/OVERNIGHT_SUMMARY.md) — overnight program through Round 8  
3. [docs/INDEX.md](docs/INDEX.md) — every doc with a one-line summary  
4. [RESEARCH_STATUS.md](RESEARCH_STATUS.md) — pre-overnight verified checks and limits  

## Reproduce checks

Requires Python 3.10+ and NumPy/SciPy.

```sh
python -m pip install -r requirements.txt
python scripts/run_checks.py
```

Scripts write results under `experiments/` (live paths for CI). Organized copies also live under `results/experiments/`.

## Layout (consolidated)

| Path | Contents |
| --- | --- |
| `docs/` | Overviews, ladder, boundary picture, INDEX |
| `docs/rounds/` | Per-rung `FFFF_*.md` notes |
| `analysis/ffz/` | FFZ Python modules + `tests/` |
| `scripts/` | Scalar/metric scans + `run_checks.py` |
| `results/` | Organized json/npz/png copies |
| `experiments/` | Pre-FFZ experiment reports + live CI results |
| `gpt_physics/` | Imported local archive (smooth_hopf, coupled_stability, …) |
| `bots/gpt/research_batches/2026-11-01/` | Earlier local work (paths kept from `main`) |
| `bots/gpt/research_batches/2026-11-21-ffz/` | Original FFZ packet checkpoints (review code moved to `analysis/ffz/`) |

## November 21 FFZ benchmark

Original packet notes: [bots/gpt/research_batches/2026-11-21-ffz/README.md](bots/gpt/research_batches/2026-11-21-ffz/README.md).  
Runnable analysis: `analysis/ffz/`. Overnight lead and ruled-out branches: [docs/OVERNIGHT_SUMMARY.md](docs/OVERNIGHT_SUMMARY.md).

Current numerical work includes free leading-N auxiliary-metric scalar checks and the interacting FFZ Newton simple-pole ledger. Passing identities do not prove physical unitarity, UV completion, experimental validity, novelty, or unification.
