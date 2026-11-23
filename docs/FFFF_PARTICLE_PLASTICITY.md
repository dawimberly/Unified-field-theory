# Particle plasticity: distortion = elastic; collapse = plastic yield

Runnable: `bots/gpt/research_batches/2026-11-21-ffz/review/ffz_particle_plasticity.py`  
Test: `.../test_ffz_particle_plasticity.py` (hooked in `scripts/run_checks.py`).  
Stacked on distortion-phase (#25).  
**No tuning.** DERIVED vs ASSUMED labeled. Includes **rate-dependent / non-Newtonian** yield.

## (1) Elastic vs plastic

From #25: \(f=(\hbar^2 n/2m)|\nabla\theta|^2+V(\theta)\).

| Response | Meaning | #25 regime |
| --- | --- | --- |
| **Elastic** | Subcritical, recoverable | distorted / incomplete slip |
| **Plastic** | Barrier escape, running phase | collapsed |

Tools (**assumed** applicability to inter-branch \(\theta\)): Josephson \(I_c\), Landau \(v_c\), Prandtl–Tomlinson, tilted washboard, Bingham / Herschel–Bulkley.

## (2) Yield threshold (robustness across \(V(\theta)\))

Tilted potential \(V_F=V_0\,v(\theta)-F\theta\). Barrier vanishes at \(F_c=V_0\max v'\) (**derived** for each shape).

| Shape | \(F_c/V_0\) |
| --- | ---: |
| Washboard \(1-\cos\theta\) | \(1\) |
| \(\cos 2\theta\) | \(1\) |
| Double well \(\tfrac14(\theta^2-1)^2\) | \(2/(3\sqrt3)\approx0.385\) |

**Robust \(O(1)\)** across shapes (**derived**). Absolute \(V_0\) **assumed** (not fixed by \(C_{\mathrm{ledger}}\)).  
Critical gradient scaling: \(|\nabla\theta|_c\sim\pi/\xi_\theta\) with \(\xi_\theta=\hbar\sqrt{n/(2m V_0)}\) (**derived** scaling).

## (3) Hysteresis (RCSJ / Stewart–McCumber)

Underdamped washboard (\(\beta_c\gg1\)): switching vs retrapping, loop width \(\approx1-4/(\pi\sqrt{\beta_c})\) (textbook approx., **assumed** prefactors).

| \(\beta_c\) | hysteretic? | loop width |
| ---: | --- | ---: |
| \(0.1\) | no | \(0\) |
| \(10\) | yes | \(\sim0.60\) |
| \(100\) | yes | \(\sim0.87\) |
| \(10^4\) | yes | \(\sim0.99\) |

**Distinct from plain DP:** DP has a single \(\tau=\hbar/E_G\), no forward/return asymmetry.

## (4) Yield length as \(R_0\)?

\(L_*=1/|\nabla\theta|_c\sim\xi_\theta/\pi\).

| Estimate | \(L_*\) | vs \(5.4\times10^{-11}\,\mathrm{m}\) |
| --- | ---: | --- |
| BEC \(\xi_{\mathrm{GP}}/\pi\) | \(\sim5.6\times10^{-8}\,\mathrm{m}\) | above bound, **system-specific** |
| Cooper \(\xi/\pi\) | \(\sim3\times10^{-8}\,\mathrm{m}\) | same |
| Tuned \(V_0\) to hit bound | possible | **tuning** (not claimed) |
| Bohr \(a_0\) | \(5.29\times10^{-11}\,\mathrm{m}\) | **coincidence only** (from #25; not derived from yield) |

**Plain: no** untuned universal \(R_0\) from yield + \(C_{\mathrm{ledger}}\).

## (5) Case signatures

| Case | \(\Omega_G=E_G/\hbar\) | Distinct tests |
| --- | ---: | --- |
| BEC \(10^6\) Rb, \(d=1\,\mathrm{mm}\) | \(\sim10^{-7}\,\mathrm{rad/s}\) (with \(\xi_{\mathrm{GP}}\)) | \(d\)-threshold; recoherence hysteresis; ramp-rate shift |
| Flux \(N_{\mathrm{eff}}=10^9\) | tiny vs electrical \(I_c\) | JJ hysteresis/KFD already exist; **gravity tilt negligible** (honest) |
| Nanosphere \(10^{-14}\,\mathrm{kg}\) | \(\sim20\,\mathrm{rad/s}\) (geometric \(R\)) | best DP-relevant; ramp separation rate vs collapse |

## Rate-dependent / non-Newtonian yield (added)

| Model | Claim | Label |
| --- | --- | --- |
| Herschel–Bulkley \(\tau=\tau_y+K\dot\gamma^n\) | \(n>1\) thickening / \(n<1\) thinning of effective yield with rate | analogy **assumed**; trend **derived** from form |
| Kurkijärvi–Fulton–Dunkleberger | Faster ramp \(\Rightarrow\) higher apparent switching threshold | **derived** washboard-escape structure; prefactors **assumed** |
| Plain DP | \(\tau=\hbar/E_G\) **rate-independent** | **derived** contrast |

**Testable signature:** collapse/recoherence threshold **shifts with preparation sweep rate** (delayed collapse at fast sweep = thickening/KFD). Plain DP does not. Absolute numbers need \(V_0,\beta_c,T\) (**assumed**). For BEC/flux, gravitational drive is tiny vs other scales (honest negative on observability there).

## Plain statement

Plasticity + hysteresis + rate-dependent yield are coherent extensions of the #25 distortion phase and give **distinct** predictions versus plain Diósi–Penrose. They still do **not** derive a Gran-Sasso \(R_0\) without tuning; the Bohr-radius near-hit stays a coincidence.
