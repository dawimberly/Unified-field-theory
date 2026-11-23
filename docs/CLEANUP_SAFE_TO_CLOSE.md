# Safe to close after this PR merges

Only after `cleanup/consolidate` is merged to `main`. This PR does **not** close PRs or delete branches.

## Contained open PRs / branches

All content from these heads is included (latest stacked tip for the newton overnight chain; plus #1 and #3). Some paths were moved (`docs/rounds/`, `analysis/ffz/`, `results/`) and links updated — that is intentional consolidation, not content loss.

| PR | Branch |
| --- | --- |
| #1 | `cursor/matter-organization-gravity-a179` |
| #2 | `cursor/ffz-newton-pole-review-2ac9` |
| #3 | `cursor/import-gpt-physics-2ac9` |
| #4 | `cursor/ffz-simple-pole-ledger-2ac9` |
| #6 | `cursor/ffz-simple-soft-linear-2ac9` |
| #7 | `cursor/ffz-rigor-a-graviton-2ac9` |
| #8 | `cursor/ffz-simple-gap-audit-2ac9` |
| #9 | `cursor/ffz-gap-physical-2ac9` |
| #10 | `cursor/ffz-b-nonoverlap-finite-2ac9` |
| #11 | `cursor/ffz-finite-eps-crosscheck-2ac9` |
| #12 | `cursor/ffz-b-literal-coupling-2ac9` |
| #13 | `cursor/ffz-dred-yukawa-2ac9` |
| #14 | `cursor/ffz-heat-kernel-a1-2ac9` |
| #15 | `cursor/ffz-superspace-ward-2ac9` |
| #16 | `cursor/ffz-condensate-r-2ac9` |
| #17 | `cursor/ffz-dynamical-condensate-2ac9` |
| #18 | `cursor/ffz-bgfield-rstar-2ac9` |
| #19 | `cursor/ffz-self-consistent-loop-2ac9` |
| #20 | `cursor/ffz-surface-term-2ac9` |
| #21 | `cursor/ffz-xi-running-2ac9` |
| #22 | `cursor/ffz-gauge-contact-2ac9` |
| #23 | `cursor/ffz-scheme-invariant-2ac9` |
| #24 | `cursor/ffz-induced-collapse-2ac9` |
| #25 | `cursor/ffz-distortion-phase-2ac9` |
| #26 | `cursor/ffz-particle-plasticity-2ac9` |
| #27 | `cursor/ffz-sm-roster-2ac9` |
| #28 | `cursor/ffz-three-phase-2ac9` |
| #29 | `cursor/ffz-two-phase-2ac9` |
| #30 | `cursor/ffz-shared-boundary-2ac9` |
| #31 | `cursor/ffz-ir-roster-2ac9` |
| #32 | `cursor/ffz-trace-anomaly-2ac9` |
| #33 | `cursor/ffz-ledger-reconcile-2ac9` |
| #34 | `cursor/ffz-lambda-pressure-2ac9` |
| #35 | `cursor/ffz-string-kk-2ac9` |
| #36 | `cursor/ffz-hadronic-a1-2ac9` |
| #37 | `cursor/ffz-ee-confinement-2ac9` |
| #38 | `cursor/ffz-running-gind-2ac9` |
| #39 | `cursor/ffz-confinement-snap-2ac9` |
| #40 | `cursor/ffz-curvature-match-2ac9` |
| #41 | `cursor/ffz-spectral-gind-2ac9` |
| #42 | `cursor/ffz-ledger-confinement-2ac9` |
| #43 | `cursor/ffz-boundary-picture-2ac9` |
| #44 | `cursor/ffz-casimir-lesage-2ac9` |
| #45 | `cursor/ffz-uv-decompose-2ac9` |
| #46 | `cursor/ffz-gut-spectrum-2ac9` |
| #47 | `cursor/ffz-round4-summary-2ac9` |
| #48 | `cursor/ffz-uv-higgsing-2ac9` |
| #49 | `cursor/ffz-planck-stress-2ac9` |
| #50 | `cursor/ffz-uv-consistency-2ac9` |
| #51 | `cursor/ffz-round5-summary-2ac9` |
| #52 | `cursor/ffz-minimal-spectrum-2ac9` |
| #53 | `cursor/ffz-positivity-region-2ac9` |
| #54 | `cursor/ffz-gap-sectors-2ac9` |
| #55 | `cursor/ffz-round6-summary-2ac9` |
| #56 | `cursor/ffz-sec-collapse-2ac9` |
| #57 | `cursor/ffz-deconf-gind-2ac9` |
| #58 | `cursor/ffz-bounce-criterion-2ac9` |
| #59 | `cursor/ffz-round7-summary-2ac9` |
| #60 | `cursor/ffz-surface-audit-2ac9` |
| #61 | `cursor/ffz-edge-modes-2ac9` |
| #62 | `cursor/ffz-phase-lock-2ac9` |
| #63 | `cursor/ffz-round8-summary-2ac9` |

(#5 was not open at consolidation time.)

## How content was taken

- Newton overnight stack: merge of tip `cursor/ffz-round8-summary-2ac9` (latest versions of all stacked files).
- #1 and #3: merged explicitly onto this branch.
- Exact duplicate files under `gpt_physics/` top-level vs nested packages removed; differing copies kept.
