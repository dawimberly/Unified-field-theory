# Penrose–Diósi collapse in induced gravity (Sakharov, \(\kappa=0\))

Runnable: `analysis/ffz/ffz_induced_collapse.py`  
Test: `.../test_ffz_induced_collapse.py` (also hooked from `scripts/run_checks.py`).  
Stacked on the scheme-invariant rung (#23).  
**No tuning.** DERIVED vs ASSUMED labeled. The Newton gap \(\approx+0.0097\) is **not** explained here — using \(C_{\mathrm{ledger}}\) vs \(1/24\) only compares two \(G_{\mathrm{ind}}\) matching choices.

## (1) Standard Diósi–Penrose

\[
\tau=\frac{\hbar}{E_G},\qquad
E_G=\hbar\,\Lambda(d),
\]

with \(\Lambda\) the collapse rate for a **uniform sphere** of mass \(M\), smearing radius \(R_0\), branch separation \(d\), \(\lambda=d/(2R_0)\):

\[
\Lambda=\frac{6GM^2}{5R_0\hbar}\times
\begin{cases}
\tfrac53\lambda^2-\tfrac54\lambda^3+\tfrac16\lambda^5 & 0\le\lambda\le1,\\
1-\dfrac{5}{12\lambda} & \lambda\ge1.
\end{cases}
\]

Far limit \(d\gg R_0\): \(\tau=5R_0\hbar/(6GM^2)\). **Derived** (textbook DP).

Proton ballpark (\(R_0\sim10^{-15}\,\mathrm{m}\), \(d\gg R_0\)): \(\tau\sim1.5\times10^7\,\mathrm{yr}\) from the exact far-limit formula (Wikipedia’s \(\sim10^6\,\mathrm{yr}\) is the same order) — **derived**.

### Gran Sasso 2020 / Donadi et al. (Nat. Phys. 2021)

The **parameter-free** DP version with \(R_0=\) nuclear positional spread (\(\sim5\times10^{-12}\,\mathrm{m}\) for Ge) is **ruled out** by the non-observation of collapse-induced radiation. Bound: \(R_0\gtrsim5.4\times10^{-11}\,\mathrm{m}\) (95% CL). **Assumed** (published experiment; not re-derived here). Honest: any DP claim that keeps nuclear-size \(R_0\) is experimentally excluded.

## (2) Induced-gravity version (\(\kappa=0\))

Packet: \(\Gamma_{\mathrm{ind}}\supset -R/(16\pi G_{\mathrm{ind}})\). The Newton simple-pole coefficient \(C\) (units \(h^2 z/K\)) matches onto that EH residue, so

\[
\frac1{G_{\mathrm{ind}}}\propto C
\quad\Rightarrow\quad
G_{\mathrm{ind}}\propto\frac1C.
\]

**Derived** from the EH matching. **Assumed:** this three-loop \(C\) dominates the Newtonian \(G\) relevant for \(E_G\) (or the comparison is only of this matching piece).

Then \(E_G\propto G\) and \(\tau\propto1/E_G\) give

\[
\frac{\tau(C)}{\tau(C_{\mathrm{ref}})}=\frac{C}{C_{\mathrm{ref}}}.
\]

| Quantity | Value | Label |
| --- | ---: | --- |
| \(C_{\mathrm{ledger}}\) | \(0.05131719\) | ledger |
| \(C=1/24\) | \(0.04166667\) | matching target |
| \(G_{\mathrm{ind}}/G_{1/24}\) | \(0.8119\) | **derived** |
| \(\tau(\mathrm{induced})/\tau(1/24)\) | \(\boldsymbol{1.2316}\) | **derived** |

Larger \(C\) ⇒ smaller \(G_{\mathrm{ind}}\) ⇒ **longer** collapse time (weaker gravity). Direction is derived, not assumed.

### Does induced gravity supply \(R_0\)?

| Scale | Size |
| --- | ---: |
| \(\ell_P(G_N)\) | \(\sim1.6\times10^{-35}\,\mathrm{m}\) |
| Nuclear size | \(\sim10^{-15}\)–\(10^{-14}\,\mathrm{m}\) (ruled out as parameter-free \(R_0\)) |
| Gran Sasso min \(R_0\) | \(\gtrsim5.4\times10^{-11}\,\mathrm{m}\) |

**Plain answer: no.** Induced UV cutoffs (\(\ell_P\), matter-loop Compton scales) do **not** naturally equal a Gran-Sasso-viable smearing length. \(R_0\) remains an extra assumption — the known weak point of DP. **Derived** scale comparison; **assumed** that DP still needs a separate \(R_0\).

Coincidence flag: none claimed between \(\ell_P\) and nuclear/\(5.4\times10^{-11}\,\mathrm{m}\) scales.

## (3) Condensate / phase-lock estimates

Induced times use \(\tau_{\mathrm{ind}}=\tau(G_N)\times C_{\mathrm{ledger}}/(1/24)\), with the **assumed** identification \(G_N\leftrightarrow\) matching at \(C=1/24\).

| Case | \(M\) | \(R_0\), \(d\) | \(\tau(G_N)\) | \(\tau(\mathrm{ind})\) | Coherence scale | Gravity dominates? |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| (a) BEC \(10^6\) \(^{87}\)Rb | \(1.44\times10^{-19}\,\mathrm{kg}\) | \(10\,\mu\mathrm{m}\), \(1\,\mathrm{mm}\) | \(6.36\times10^{8}\,\mathrm{s}\) | \(7.84\times10^{8}\,\mathrm{s}\) | \(\sim1\,\mathrm{s}\) (assumed) | **No** |
| (b) Flux qubit / Cooper pairs \(N_{\mathrm{eff}}=10^9\) | \(1.82\times10^{-21}\,\mathrm{kg}\) | \(1\,\mu\mathrm{m}\), \(1\,\mu\mathrm{m}\) | \(1.49\times10^{12}\,\mathrm{s}\) | \(1.84\times10^{12}\,\mathrm{s}\) | \(\sim50\,\mu\mathrm{s}\) (assumed) | **No** |
| (c) Nanosphere \(10^{-14}\,\mathrm{kg}\) | \(10^{-14}\,\mathrm{kg}\) | \(R_0=1.06\,\mu\mathrm{m}\) from \(\rho=2\times10^3\), \(d=R_0\) | \(5.3\times10^{-2}\,\mathrm{s}\) | \(6.5\times10^{-2}\,\mathrm{s}\) | \(\sim1\,\mathrm{ms}\) (assumed) | **No** at \(d=R_0\) with that coherence; closer if coherence longer / \(d\) larger |

Labels: \(\tau\) from DP **derived**; \(R_0\), coherence scales, \(N_{\mathrm{eff}}\), \(\rho\) **assumed**. Real flux-qubit mass displacement may be smaller ⇒ even larger \(\tau\).

## (4) Phase framing

\[
\frac{E_G\,t}{\hbar}\sim1
\quad\Longleftrightarrow\quad
t\sim\tau=\frac{\hbar}{E_G}.
\]

This is a **rewrite** of the same DP criterion in phase/frequency language (relative gravitational phase between branches). It does **not** add new dynamics by itself; a concrete stochastic law (Diósi Lindblad) is an extra **assumed** step. **Derived** equivalence of the two writings.

## Relation to the \(+0.0097\) Newton gap

None derived. Using \(C_{\mathrm{ledger}}\) rather than \(1/24\) only rescales \(\tau\) by \(\approx1.232\). That does not explain why the ledger differs from \(1/24\).

## Plain statement

In the \(\kappa=0\) induced-gravity reading, collapse is slightly **slower** than with a \(1/24\)-matched \(G\) (\(\tau\) ratio \(\approx1.23\)). Induced gravity does **not** fix DP’s \(R_0\) problem (Gran Sasso still applies). Macroscopic BEC/flux estimates stay far from gravity-dominated collapse; a \(10^{-14}\,\mathrm{kg}\) nanosphere is the only ballpark where DP times can approach coherence scales, and only with assumed \(R_0\) and isolation.
