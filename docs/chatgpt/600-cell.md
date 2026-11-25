# 600-cell geometry

**Conversation date range** 2026-10-13 to 2026-10-25 (archive dates; original transcript dates UNVERIFIED)

**ASSUMED** The 600-cell was considered as a finite network for studying phase locking, topology and information flow. **CITED** Coxeter, *Regular Polytopes*, supplies the standard geometric counts summarized below; they are not a new derivation. **ASSUMED** The proposed physical use places states on vertices and interactions on edges. **UNVERIFIED** No emergence of spacetime or physical force law was established.

## Geometry

**CITED** Coxeter, *Regular Polytopes*: The 600-cell, or hexacosichoron, is a regular four-dimensional polytope with Schläfli symbol $\{3,3,5\}$. Its boundary triangulates $S^{3}$, rather than $S^{2}$.

| Component | Count |
|---|---:|
| **CITED** Vertices (Coxeter, *Regular Polytopes*) | 120 |
| **CITED** Edges (Coxeter, *Regular Polytopes*) | 720 |
| **CITED** Triangular faces (Coxeter, *Regular Polytopes*) | 1,200 |
| **CITED** Tetrahedral cells (Coxeter, *Regular Polytopes*) | 600 |

**UNVERIFIED** Each vertex has 12 neighbors and an icosahedral link. On the unit sphere, neighboring vertices have dot product $\varphi/2\approx0.809$, where $\varphi=(1+\sqrt5)/2$. The archived coordinate construction uses eight axis points, sixteen points with every coordinate equal to $\pm1/2$, and ninety-six even permutations of $(0,\pm1/2,\pm\varphi/2,\pm1/(2\varphi))$. The coordinate set was identified with the binary icosahedral group; no independent verification is recorded here.

## Proposed physics use

**ASSUMED** Vertices would carry phases, quantum states or other field variables; edges would carry interactions. The finite, symmetric network could test whether local updates synchronize, drift or mix. Selecting this network supplies a geometry at the outset; it does not derive geometry from microscopic physics or establish that physical spacetime has this lattice.

## Separate archived market application

**UNVERIFIED** The archived discussion describes `hypersphere_600cell_backtest.py` using the 120 vertices as weight directions over four lagged features: stock momentum, market permutation entropy, signed sector synchronization and signed sector-to-stock transfer entropy. It reports a 720-edge graph check, but the replay scored feature vectors by a dot product rather than evolving states along edges. Nonpositive scores stayed in cash.

**UNVERIFIED** The reported first-window selection used 75% of cached history. Cell 98 had approximate weights $(0.5,0,0.31,0.81)$, with reported returns of $+41.46\%$ and $+50.17\%$ on the selection and later windows, compared with an MA50 baseline of $-6.85\%$ and $+42.50\%$. These figures have not been reproduced in this archive and provide no physics evidence. The holdout had already been inspected and was to remain locked.

## Limits

**UNVERIFIED** A reproducible physics calculation would need an explicit state space, action or update law, observables and refinement checks. No such physics calculation is supplied in this topic; the standard geometric counts are sourced to Coxeter, *Regular Polytopes*. The market replay and the proposed physics network are separate applications.
