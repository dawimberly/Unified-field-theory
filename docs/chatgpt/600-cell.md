# 600-cell geometry

**Conversation date range** 2026-10-25 (archive discussion; original transcript date UNVERIFIED)

**ASSUMED** The 600-cell was considered as a geometric candidate for organizing possible motion. **UNVERIFIED** The discussion records exploratory correspondences. **UNVERIFIED** A physical derivation from the geometry remains open.

## Research record

# Cursor Conversation Export

---

## 1. USER

**UNVERIFIED** <timestamp>Tuesday, Oct 13, 2026, 8:24 PM (UTC-5)</timestamp>
**UNVERIFIED** You're thinking of the 600-cell, also called the hexacosichoron.

**UNVERIFIED** [File:Schlegel wireframe 600-cell vertex-centered.png - Wikimedia Commons](https://images.openai.com/static-rsc-4/z4d2EI8wlscJJGGBX6sjBIb8bfZ6qrM2Q6liM7gDug0sFTx_vQIH659aNnbJdSLKmsnMuOsn9yfFhA7b_ptE8hYhfVgX7JBab-NsA94rEUIBs47YcHGqHa9IRwkyDQjF3Tx7RNAjFtJca72B9-2Vv5VUksrqh0MQRoD9BmhQG14?purpose=inline)

**UNVERIFIED** [File:600-cell t02 B3.svg - Wikimedia Commons](https://images.openai.com/static-rsc-4/aVwBm8r8Y3zSlfcBsiZgUbwX4ndQYd6P7hrAHeV2ANXV846dLzgZivlqa8YuVp88_BQZimVGE9qaSb-VMWTjJ1Nh0XpAICX94N4h_Iv9XdWJXh1O2YvCr1B-2lDTLXINYBh_v68JT9GPFZZE8SkTDzsP-vy3j8gdKVH3wOM3KUE?purpose=inline)

**UNVERIFIED** 7

**UNVERIFIED** It's a regular four-dimensional polytope whose boundary forms a triangulation of the three-dimensional sphere, \$S^3\$.

**UNVERIFIED** Here are its properties:

| Component           | Number |
| ------------------- | ------ |
| **DERIVED** Vertices (nodes)    | 120    |
| **DERIVED** Edges (connections) | 720    |
| **DERIVED** Triangular faces    | 1,200  |
| **DERIVED** Tetrahedral cells   | 600    |

### Why we were interested in it

**UNVERIFIED** We discussed using the 600-cell as a discrete starting grid for your theoretical physics framework.

**UNVERIFIED** Instead of assuming that space is fundamentally continuous, we could explore a network where:

- **UNVERIFIED** Vertices carry quantum states, phases, or information.
- **UNVERIFIED** Edges represent relationships and interactions.
- **UNVERIFIED** Tetrahedral cells provide a discrete geometric structure.
- **UNVERIFIED** Collective behavior could potentially generate effective geometry.

**UNVERIFIED** The 600-cell is mathematically well-defined, but there is no established evidence that spacetime is physically constructed from this particular lattice.

**UNVERIFIED** One important detail: the 600-cell is not an ordinary sphere with 600 points. It has 120 vertices and 600 tetrahedral cells, and its boundary is topologically \$S^3\$, rather than the familiar two-dimensional spherical surface \$S^2\$.

**UNVERIFIED** We originally considered it because it could provide a manageable mathematical laboratory for studying phase locking, topology, and emergent geometry.

**UNVERIFIED** And given our recent discussion about entropy and the stock market, it might also make an interesting test network for information propagation and synchronization.&#x20;

---

## 2. ASSISTANT

**UNVERIFIED** The 600-cell is a 120-vertex, 600-tetrahedron triangulation of $S^3$, so it can be a finite network for phase locking and information flow. I’ll check how that lines up with the entropy and market work already in this project.

---

## 3. ASSISTANT

**UNVERIFIED** The counts are right: the 600-cell is the regular 4-polytope with Schläfli symbol {3,3,5}. It has **120 vertices**, **720 edges**, **1,200 triangles**, and **600 tetrahedra**, and that boundary is a triangulation of $S^3$.

**UNVERIFIED** Each vertex meets 12 edges, and the link of a vertex is an icosahedron. With the vertices placed on the unit 3-sphere, nearest neighbors sit at dot product $\varphi/2 \approx 0.809$, where $\varphi$ is the golden ratio. The coordinate set is the binary icosahedral group: the 8 axis points, the 16 points with coordinates $\pm 1/2$, and the 96 even permutations of $(0,\ \pm 1/2,\ \pm \varphi/2,\ \pm 1/(2\varphi))$.

**UNVERIFIED** That object is already in the research replay, used in one specific way. `hypersphere_600cell_backtest.py` builds those 120 unit vectors and checks that the nearest-neighbor graph has 720 edges. Each vertex is a weight direction over four lagged features:

- **UNVERIFIED** **M** — stock momentum at the 10:15 bar
- **UNVERIFIED** **O** — market permutation entropy through the prior close
- **UNVERIFIED** **C** — signed sector synchronization
- **UNVERIFIED** **F** — signed sector-to-stock transfer entropy

**UNVERIFIED** The score is the dot product. A non-positive score stays in cash. On the first 75% of the cached history, the leading directions put most of their weight on momentum and information flow. Cell 98, weights about $[0.5,\ 0,\ 0.31,\ 0.81]$, returned +41.46% in that window and +50.17% on the later holdout, against an MA50 baseline of −6.85% then +42.50%. That holdout has already been read, so it should stay locked.

**UNVERIFIED** That replay treats the vertices as directions in the four-feature space. It does not put a state on each vertex and pass information along the 720 edges. A propagation test would use the 12-regular graph itself: 120 nodes, each with a phase or a market state, updates only along the icosahedral neighborhood, and a readout for whether the closed $S^3$ topology locks, drifts, or mixes. The symmetry makes that a clean laboratory. It is still a model of the network, not a claim that spacetime or the market is this lattice.

---
