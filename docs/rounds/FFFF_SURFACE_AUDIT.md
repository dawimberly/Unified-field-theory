# Boundary / surface / total-derivative audit (AG)

Runnable: `analysis/ffz/ffz_surface_audit.py` + test. Stacked on #59.  
**No tuning.** Fields-only. Attack GAP \(+0.009651\) with terms internal to original FFZ fields.

## Setup

Three-loop WZ+Yukawa Newton residue in units \(h^2 z/K\): weakly curved \(\mathbb{R}^4\), compactly supported metric probe, **no geometric \(\partial M\)**.

## Results

| Term | \(\delta C\) | Scheme-independent into \(h^2z/K\)? | Label |
| --- | --- | --- | --- |
| \(\Box R\) | **0** | yes (vanishes) | **derived** (Stokes) |
| GHY / extrinsic \(K\) | **0** | yes (no \(\partial M\)) | **derived** |
| Euler \(E_4\) | **0** | topology SI; **not** in \(h^2z/K\) | **derived** (Lovelock) |
| \(a_2\)-type TD / improvements | **0** | yes (wrong HK grade / channel) | **derived** |
| Finite vacuum EH \(\kappa\) | \(\kappa=-\mathrm{GAP}\) if chosen | **no** (scheme) | allowed CT (#23) |

Derived surface/TD sum \(=0\). GAP stays \(+0.009651\pm0.0024\).

## Plain

Nothing was “dropped” that would scheme-independently shift the simple pole. Euler is topological / does not vary to Einstein. Absorbing \(\kappa=-\mathrm{GAP}\) is the same scheme contact already named in #23 — not a missing bulk topology.
