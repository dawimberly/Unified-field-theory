# Quantum-to-Einstein gravity research

An exploratory research notebook, not a demonstrated unified theory. Start with [RESEARCH_STATUS.md](RESEARCH_STATUS.md) for current results, limits, and next calculations.

## November 21 FFZ benchmark

[Review packet for the interacting Newton matching calculation](bots/gpt/research_batches/2026-11-21-ffz/README.md). It includes the current pole checks and the unresolved FFff graphs and counterterms.

## Reproduce the current checks

Requires Python 3.10+ and NumPy/SciPy.

```sh
python -m pip install -r requirements.txt
python scripts/run_checks.py
```

Scripts write results under experiments/ relative to the repository, independent of the working directory. Existing reports describe the original exploratory runs; scripts are the executable reproduction source. See each report for conventions and scope.

## Contents

- `scripts/`: full scalar-response scan, source projection, metric tensor/contact Ward check, continuum spectral check.
- `experiments/`: calculation reports and numerical results.
- [RESEARCH_STATUS.md](RESEARCH_STATUS.md): verified checks, assumptions, unresolved issues, next work.
- [docs/POSSIBILITY_CLOUD_ORIGINAL.md](docs/POSSIBILITY_CLOUD_ORIGINAL.md): earlier qualitative Possibility Cloud proposal, preserved verbatim. Its claims are hypotheses, not established conclusions of the current calculations.
- [docs/GRAVITON_ROADMAP.md](docs/GRAVITON_ROADMAP.md): how an induced \(R\) term relates to a graviton kinetic term, the \(1/24\) matching, three graviton conditions, and what remains open.

Current numerical work concerns a free, leading-N, auxiliary-metric scalar determinant. Equivalence to constrained CEV is not established here. Passing numerical identities does not prove physical unitarity, UV completion, experimental validity, novelty, or unification.

## November 1 field and dimension follow-up

[Latest discussion, reports, and computation archive links](bots/gpt/research_batches/2026-11-01/README.md). The field model is a separate exploratory line from the CEV calculations. Its best tested configuration is not a certified full equilibrium. Raw-file transfer coverage is explicitly documented in the batch.
