# Proven-methods ladder for the FFZ +0.01 simple-pole gap

Dan’s plan: try established peer-reviewed methods in turn; if none closes the
\(\sim+0.01\) projected simple-pole gap vs \(1/24\), a new idea is needed.

Prescription throughout: \(\xi=1/6\), \(a_X=0\), \(a_Y=-1\), \(L=0\), units \(h^2 z/K\).

| # | Method | Status | Result (concrete) | PR |
| ---: | --- | --- | --- | --- |
| 1 | DREG region / projective diagram ledger (FFff overlap fix \(\Delta a_1=-16\), \(\Delta c_0=64/3\)) | **done** | Higher poles cancel; simple total \(\approx+0.0514\); **gap \(\approx+0.0097\pm0.002\)** | [#4](https://github.com/dawimberly/Unified-field-theory/pull/4)–[#10](https://github.com/dawimberly/Unified-field-theory/pull/10) |
| 2 | Independent finite-\(\varepsilon\) Laurent fit (no stored \(c_0/a_1\)) | **done** | Gap **\(+0.00954\pm0.0036\)** survives | [#11](https://github.com/dawimberly/Unified-field-theory/pull/11) |
| 3 | DRED + global Yukawa weight + BPHZ/R | **done** | Higher poles cancel only after R (soft-soft \(-16,64/3\)); **gap \(+0.00975\pm0.0024\)** stands | [#13](https://github.com/dawimberly/Unified-field-theory/pull/13) |
| 4 | Heat-kernel / Seeley–DeWitt \(a_1\) route to induced \(R\) (BG effective action) | **partial** | Full 3-loop interacting HK **infeasible** here; lower-loop HK inputs that fix \(c_Y=1/24\) and CT **verified**; gap **\(+0.00975\)** unchanged | [#14](https://github.com/dawimberly/Unified-field-theory/pull/14) |
| 5 | Superspace / supergraph method | **partial** | SUSY forces \(c_Y=1/24\), **not** baseline\(=1/24\); Ward inventory PASS; no missing partner; gap **\(+0.00975\)** stands | [#15](https://github.com/dawimberly/Unified-field-theory/pull/15) |
| 6 | Background-field method with \(R^*\) operation | **partial** | R* IR ≡ soft forest IE; soft-linear not \(1/\varepsilon\); no missed IR CT; total \(+0.05132\pm0.0024\); **gap \(+0.00965\pm0.0024\)**; pieces unchanged | [#18](https://github.com/dawimberly/Unified-field-theory/pull/18) |
| ★ | **New idea:** Higgs-like condensate / phase locking (\(\xi v^2\) \(R\) contact) | **done (negative)** | Soft/F-term vev **assumed**; tree \(\delta=K v^2/(6 h^2 z)\); gap-absorbing \(v^2/z=3.70\times10^{-4}\,h^2\) **not predicted**; free parameter / tuning | [#16](https://github.com/dawimberly/Unified-field-theory/pull/16) |
| ★★ | **New idea:** Dynamical condensate (BCS/NJL / soft gap equation) | **done (negative)** | SUSY forbids \(v\neq0\); soft/NJL needs free \(m_s^2/z\) or \(\Lambda/\sqrt{z}\); ratio pred/target \(=0\) or \(O(10^{-2}\)–\(10^6)\) — not 1 without tuning | [#17](https://github.com/dawimberly/Unified-field-theory/pull/17) |
| ★★★ | **New idea:** Self-consistent loop \(v\leftrightarrow G_{\mathrm{ind}}\leftrightarrow m_s\) (gravity mediation) | **done (negative)** | \(F\sim z\): nontrivial fixed point but \(\delta_R/\mathrm{GAP}\sim87\) (\(h=1\)); \(F=\|y\|v^2\): only \(v=0\); free \(\varphi\) restores tuning | *(this branch → stacked PR)* |

## How to read status

- **done** — method executed; concrete numbers reported; gap still open unless stated otherwise.
- **partial** — only the feasible subset of the method was run; the blocked part is named in the method’s doc.
- **next / queued** — not started.
- **new idea** — not an established peer-reviewed method; exploratory. **done (negative)** = executed and does not close the gap without tuning.

## Gap summary (standing)

\[
\text{ledger total}\approx +0.0514,\qquad
\frac1{24}\approx +0.041667,\qquad
\text{GAP}\approx +0.0097\pm 0.002.
\]

**Established ladder (1–6) exhausted without closing the gap.** Condensate new-ideas also negative. A genuinely new approach is needed.

## Pointers

- Ledger / audits: `bots/gpt/research_batches/2026-11-21-ffz/review/`
- Heat-kernel rung: `FFFF_HEAT_KERNEL_A1.md`, `ffz_heat_kernel_a1.py`
- Superspace / Ward rung: `FFFF_SUPERSPACE_WARD.md`, `ffz_superspace_ward.py`
- BG-field / R*: `FFFF_BGFIELD_RSTAR.md`, `ffz_bgfield_rstar.py`
- Condensate / new idea: `FFFF_CONDENSATE_R.md`, `ffz_condensate_r.py`
- Dynamical condensate: `FFFF_DYNAMICAL_CONDENSATE.md`, `ffz_dynamical_condensate.py`
- Self-consistent loop: `FFFF_SELF_CONSISTENT_LOOP.md`, `ffz_self_consistent_loop.py`
- Graviton context: `docs/GRAVITON_ROADMAP.md`
