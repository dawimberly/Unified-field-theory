# Local GPT physics archive (import)

Imported from Dan’s local folder
`C:\Users\Owner\OneDrive\Documents\ChatGPT\Physics`
(paths normalized from Windows backslashes). Snapshot date context: late Oct–early Nov 2026 write-ups around Candidate 0, Hopf textures, and gravity cross-checks.

Some scripts and markdown also appear both at the top level and inside the packaged folders; both copies are kept as in the source tree.

## Write-up index

| Write-up | Covers | Scripts / data |
| --- | --- | --- |
| [`gravity_cross_reference_2026_10_30.md`](gravity_cross_reference_2026_10_30.md) | Cross-reference of gravity / Higgs / Hopf / CEV threads; no demonstrated unification. Coverage caveats and Higgs-sector consistency test. | [`gravity_cross_reference_checks.py`](gravity_cross_reference_checks.py), [`gravity_cross_reference_results.json`](gravity_cross_reference_results.json), [`gravity_review_coverage.json`](gravity_review_coverage.json), [`scalar_rational_certificate_rerun.json`](scalar_rational_certificate_rerun.json), [`tensor_rouche_certificate_rerun.json`](tensor_rouche_certificate_rerun.json). Extra prior certificates and Hopf reference inputs only in [`gravity_cross_reference_2026_10_30.zip`](gravity_cross_reference_2026_10_30.zip) (`prior_certificates/`, `reference_inputs/`). |
| [`field_transition_cross_reference_2026_11_01.md`](field_transition_cross_reference_2026_11_01.md) | Assessment of a “crest / field-transition” collapse idea against existing Candidate 0 and project threads; proposes a two-location single-particle test. | Conceptual note; no dedicated companion script in this import. |
| [`dimension_gravity_investigation_2026_11_01.md`](dimension_gravity_investigation_2026_11_01.md) | Mechanism audit of why \(d=3\) is special for fluids/gravity vs project scaling; dilation balance for \(2<d<4\). | [`dimension_mechanism_checks_2026_11_01.py`](dimension_mechanism_checks_2026_11_01.py), [`dimension_mechanism_checks_2026_11_01.json`](dimension_mechanism_checks_2026_11_01.json). |
| [`coupled_stability_2026_11_01/REPORT.md`](coupled_stability_2026_11_01/REPORT.md) | Coupled \(\rho\)+\(n\) lattice relaxation: does **not** establish a stable localized continuum knot; winding diagnostics and vacuum spectrum. | Package folder [`coupled_stability_2026_11_01/`](coupled_stability_2026_11_01/): `coupled_stability_test_*.py`, `coupled_hessian_check_*.py`, `coupled_vacuum_spectrum_*.py`, `knotted_cloud_dynamics.py`, `package_coupled_stability_*.py`; JSON (`N32/N48_results`, Hessian/force/vacuum checks, `manifest.json`); `.npz` candidates; `relaxation_energy.png`. Top-level copies of several of the same `.py` files also exist. Candidate 0 action text: [`coupled_stability_2026_11_01/CANDIDATE_ACTION.md`](coupled_stability_2026_11_01/CANDIDATE_ACTION.md). |
| [`smooth_hopf_followup_2026_11_01.md`](smooth_hopf_followup_2026_11_01.md) (also under [`smooth_hopf_2026_11_01/`](smooth_hopf_2026_11_01/smooth_hopf_followup_2026_11_01.md)) | Smooth, topology-preserving Hopf trial family on \(\mathbb{R}^3\); variational energy drop vs radial control; omitted-mode / derivative checks. | Package folder [`smooth_hopf_2026_11_01/`](smooth_hopf_2026_11_01/): `smooth_hopf_variational_*.py`, `smooth_hopf_extended_*.py`, derivative / omitted / radial / smallest-mode checks, `package_smooth_hopf_*.py`; JSON summaries and order-parameter dumps; small `.npz` parameter files; `family_energy.png`. Top-level copies of the main scripts and markdown also exist. Candidate 0 action text: [`smooth_hopf_2026_11_01/CANDIDATE_ACTION.md`](smooth_hopf_2026_11_01/CANDIDATE_ACTION.md). |

## Omitted from this import

- **`coupled_stability_2026_11_01.zip`** — duplicate of the unpacked `coupled_stability_2026_11_01/` tree.
- **`smooth_hopf_followup_2026_11_01.zip`** — byte-identical to the unpacked `smooth_hopf_2026_11_01/` tree (28/28 files matched).
- **`.npz` files** — none omitted; all seven are under ~5 MB (largest: `coupled_stability_2026_11_01/N48_stationary_candidate.npz` ≈ 3.1 MB).

Kept zip with content not otherwise present: **`gravity_cross_reference_2026_10_30.zip`**.
