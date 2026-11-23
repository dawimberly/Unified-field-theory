# Graviton roadmap (status and open questions)

Plain checkpoint of how the present induced-gravity calculation relates to a graviton. Nothing here claims a completed unified field theory.

## Where we are

The active benchmark is a Wess–Zumino / chiral matter theory on an **external** metric, fixed \(\xi=1/6\), with a three-loop Newton matching test in \(d=4-2\varepsilon\). The improvement-independent simple pole of the curvature response is required to equal \(1/24\) (units \(h^2 z/K\)) at the projected point \(a_X=0\), \(a_Y=-1\).

**Status of that test (honest):**

- Triple and double poles cancel after the FFff(B) two-log overlap assembly fix (\(\Delta a_1=-16\), B double \(-1/288\)).
- That fix’s soft-soft \(W/U\) finite \(\Delta c_0=64/3\) is derived from the same overlap density; soft-linear is excluded by sector inclusion-exclusion (see review `FFFF_B_OVERLAP_RIGOR.md`).
- The projected simple pole versus \(1/24\) is **still open** (gap \(\approx+0.0097\)). Physical audit: \(c_Y=1/24\) is RG-fixed, but baseline\(=1/24\) is a matching test — residual allowed / scheme-absorbable; no missing Wick topology; see `FFFF_GAP_PHYSICAL.md`.

Passing \(1/24\) would fix the coefficient of an induced Einstein–Hilbert term from this matter sector. It would **not** by itself produce a dynamical graviton.

## From induced \(R\) to a graviton kinetic term

Write the induced low-energy piece schematically as

\[
\Gamma_{\text{ind}}[g]
\supset
-\frac{1}{16\pi G_{\text{ind}}}\int\mathrm{d}^4x\,\sqrt{g}\,R
+\cdots
\]

Expand \(g_{\mu\nu}=\eta_{\mu\nu}+\kappa h_{\mu\nu}\) with \(\kappa\sim\sqrt{G_{\text{ind}}}\). At quadratic order in \(h\) one obtains the linearized Einstein kinetic operator for \(h_{\mu\nu}\) (gauge-fixed or in terms of gauge-invariant projectors). That is the sense in which an induced \(R\) term “yields a graviton”: the metric fluctuation propagates with the Einstein kinetic term whose residue is set by \(G_{\text{ind}}\).

If the matched three-loop simple pole equals \(1/24\), that number enters the matching condition that determines \(G_{\text{ind}}\) relative to the matter couplings and mass scale in this scheme (together with lower-loop and threshold data already in the packet). A failed \(1/24\) test means this particular matter-induced coefficient does not reproduce the required Newton residue; it does not by itself invent a different spin-2 particle.

## Three conditions any graviton candidate must pass

1. **Massless spin-2 with exactly two helicities; no ghosts / extra modes.**  
   The quadratic kinetic operator must project onto the massless spin-2 representation only. Extra scalar or vector metric modes with wrong-sign or light residues are fatal.

2. **Universal coupling / equivalence principle.**  
   The residue that defines \(G_{\text{ind}}\) must couple to the conserved stress tensor in a species-independent way at leading order (same \(G\) for all probes, up to controlled higher-dimension corrections).

3. **Weinberg–Witten evasion for composites.**  
   If the graviton is composite (emergent from matter), it must evade the Weinberg–Witten obstruction (no Lorentz-covariant conserved stress tensor for a massless spin-2 charge in the usual flat-space sense). Typical evasions: the metric is not a Poincaré-covariant local composite in flat space, or diffeomorphism invariance is built in from the start so the “composite” is the geometry itself.

The present calculation tests a necessary coefficient in (the matter-side of) (1)–(2) on an external metric. It does not yet address (3), nor does it promote \(g_{\mu\nu}\) to a path-integral variable.

## Open question: what makes the metric dynamical?

The FFZ benchmark treats \(g_{\mu\nu}\) as an external probe. A unified theory needs microscopic variables whose effective description includes a dynamical metric (or an equivalent spin-2 mode) with the properties above. Candidate directions in this repo’s history—auxiliary-metric large-\(N\), constrained measure / CEV, composite condensates—remain **open** and are not established by the three-loop pole work.

**Next calculational gate:** continue the proven-methods ladder (`docs/PROVEN_METHODS_LADDER.md`) — superspace/Ward feasible subset done (SUSY forces \(c_Y\) not baseline; gap stands); next is background-field \(R^*\) — until the \(\sim+0.01\) gap closes or every established method is exhausted. **Next theoretical gate after that:** specify the microscopic variables that make \(g_{\mu\nu}\) dynamical and check the three conditions above in that setting.
