# Results

Organized copies of numerical artifacts (`.json`, `.npz`, `.png`, …).

| Path | Source |
| --- | --- |
| `results/experiments/` | Copies from `experiments/*` (CI still writes/reads under `experiments/`) |
| `results/gpt_physics/` | Artifacts from the `gpt_physics/` import |
| `results/ffz/` | Any leftover FFZ review artifacts |
| `results/CONSOLIDATION_MANIFEST.json` | File moves performed in the cleanup PR |

Do not delete `experiments/*/results.json` — reproduction scripts depend on those live paths.
