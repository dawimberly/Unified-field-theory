# Independent screening audit

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records independent screening audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Six original isolated-sphere runs reproduce the reported values. An independent piecewise adaptive SciPy collocation solver confirms the approximation error and the effect of radial mass distribution. This is a numerical test of an $n=1$ chameleon scalar model, not an experimental detection, quantum-gravity calculation, or unification result.

## Independent results

| mu | Uniform exterior charge | Layered exterior charge | Layering reduction | Thin-shell excess (uniform / layered) |
|---|---:|---:|---:|---:|
| **DERIVED** 0.1 | 0.46179560 | 0.42491187 | 7.99% | 25.82% / 36.74% |
| **DERIVED** 0.2 | 0.11968499 | 0.11326649 | 5.36% | 21.37% / 28.24% |
| **DERIVED** 0.3 | 0.05168994 | 0.04983739 | 3.58% | 24.90% / 29.54% |

## Verification

- **DERIVED** Original versus independent exterior-charge discrepancy: at most 0.0817%.
- **DERIVED** Independent integrated exterior flux balance: maximum normalized discrepancy $1.85\times10^{-09}$.
- **DERIVED** All 62 recorded independent checks pass (solver convergence, positivity, equal continuous mass, flux balance, and refinement comparisons).
- **DERIVED** Tighter solver tolerance, doubled exterior radius, and a linear asymptotic Robin boundary were tested for both $\mu =0.1$ profiles. Each matches its baseline charge within 0.02%.
- **DERIVED** Additional original-grid and boundary checks are in extra_checks.json.

| mu | Core fraction | Intervals | Outer radius | Exterior charge |
|---|---|---:|---:|---:|
| **DERIVED** 0.1 | 0 | 192000 | 120 | 0.461984200 |
| **DERIVED** 0.1 | 0.8 | 192000 | 120 | 0.424951667 |
| **DERIVED** 0.3 | 0.8 | 192000 | 120 | 0.049856497 |
| **DERIVED** 0.1 | 0.8 | 192000 | 240 | 0.424992071 |

## Interpretation and issues

**DERIVED** The local surface flux is not the exterior Yukawa amplitude. The exterior equation contains a distributed scalar source: $[r^2u\prime]_1^r=\int_1^r ks^2(1-u^{-2})\,ds$. Independently integrating this source verifies the flux change. The mismatch is expected physics within this model, not by itself a coding error.

**DERIVED** The layered body has core radius 0.8, core density 1500, mantle density 475.409836, and mean density 1000. Uniform and layered cases have the same continuous mass and radius. This probes radial density arrangement, not chemical composition, piezoelectricity, or intrinsic spin.

**DERIVED** The original node-assigned density jumps produce a finite-grid interface error. Independent domain matching treats interfaces exactly; refined original runs approach its charge. Surface derivatives are more sensitive to this error than far-field charges.

**DERIVED** Inspection found a metadata bug: finiteWallAndExterior.exteriorBoundaryRadius is hardcoded to 40 even when the actual numerical boundary is 120. The solver uses the requested radius; the JSON label is wrong. The original core-interface alignment check also uses coreFraction*N instead of coreFraction/h. Neither invalidates these aligned isolated-sphere runs.

**DERIVED** This audit independently verifies the isolated-body calculation only. Original runs also execute chamber/wall cases, but those geometries and the other cosmology/MICROSCOPE scripts have not received the same independent verification here. No observational exclusion follows from these tests.

## Next calculation and questions

**ASSUMED** The most useful next calculation is a controlled density-contrast/core-radius scan at fixed total mass and radius, followed by matching any proposed material mechanism to an explicit scalar coupling. A scalar screening effect must be distinguished from a change in the Einstein spin-2 gravitational coupling.

1. **DERIVED** Is the target a density-arrangement dependence of an extra scalar force, or a composition/spin dependence at identical density? The latter requires new variables and a stated coupling.
2. **DERIVED** Is there a particular material pair or apparatus to model? Without one, the next scan should remain dimensionless rather than be presented as an experiment prediction.

## Reproduction

**DERIVED** Requires Python with NumPy/SciPy and Node.js. Run sequentially to avoid concurrent snapshot writes.

```powershell
python C:\workspace\physics\scripts\screening_independent_audit.py --originals
python C:\workspace\physics\scripts\screening_independent_audit.py
python C:\workspace\physics\scripts\screening_finish_audit.py
```

**DERIVED** Input scripts and the report are archived unchanged in inputs/; hashes are in input_manifest.json. Results and independent scripts are included. Original files were not edited.
