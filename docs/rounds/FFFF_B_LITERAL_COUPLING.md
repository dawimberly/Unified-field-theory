# FFff(B) literal coupling-weight cross-check

Runnable: `ffz_ffff_B_literal_coupling.py`.  
**No forest/region inclusion-exclusion. No tuning to \(1/24\) or \(-1/288\).**

Independent check of the overlap fix \(\Delta a_1=-16\), \(\Delta c_0=64/3\), which every prior method assumed.

## (b) Coupling dimension from the bare Lagrangian

In \(d=4-2\varepsilon\):

| Quantity | Value |
| --- | --- |
| \([\phi]=(d-2)/2\), \([\psi]=(d-1)/2\) | canonical |
| Yukawa \(y\,\phi\,\psi\psi\) \(\Rightarrow [y]\) | \(\varepsilon\) |
| Quartic from auxiliaries \([h]=[|y|^2]\) | \(2\varepsilon\) |
| Field variation weight \(\lambda_F=(4-d)/2\) | \(\varepsilon\) |

Local Weyl source: \(\delta y=-\varepsilon\sigma(x_v)y\) with \(\sigma\sim x_v^2\) inserts \(-\varepsilon\sum T_{ii}\) at Yukawa vertices. The fermion kinetic / heat-kernel field-weight vertex piece is \(-\lambda_F\sum T_{ii}\). Therefore as operators at this order

\[
A_{\text{coupling}}\equiv A_{\text{vertex}}\qquad(\varepsilon=\lambda_F).
\]

That is an identity, not a fit. A *literal* Lagrangian that keeps both the kinetic Weyl variation and the local Yukawa weight therefore adds \(A_{\text{vertex}}\) **twice globally**, unless a renormalized-field convention absorbs one copy into the other.

## (a) Sector residues and graph poles

Vertex-piece soft-soft residues for three weight choices (same operator, different counting):

| Weight | \(w_f\) | \(w_c\) | ll | pl | \(W/U\) \(c_0\) |
| --- | ---: | ---: | ---: | ---: | ---: |
| field_only | 1 | 0 | \(-16\) | \(-48\) | \(+64/3\) |
| coupling_only | 0 | 1 | \(-16\) | \(-48\) | \(+64/3\) |
| literal global | 1 | 1 | \(-32\) | \(-96\) | \(+128/3\) |

Graph poles after applying each assembly to the field-only projective base:

| Case | triple | double | simple | proj. double residual |
| --- | ---: | ---: | ---: | ---: |
| field_only | \(+0.020833\) | \(-0.010426\) | \(-0.052302\) | \(-6.95\times10^{-3}\) |
| literal_global (\(\Delta a_1=-64\), \(c_0{+}64/3\) only) | \(+0.020833\) | \(+0.017352\) | \(+0.054179\) | \(+2.08\times10^{-2}\) |
| literal_global (+WU + soft-linear) | \(+0.020833\) | \(+0.017352\) | \(+0.099884\) | \(+2.08\times10^{-2}\) |
| **overlap_fixed** (\(\Delta a_1=-16\), \(\Delta c_0=64/3\)) | \(+0.020833\) | \(\boldsymbol{-0.003481}\) | \(\boldsymbol{-0.032626}\) | \(\boldsymbol{-8.9\times10^{-6}}\) |

Notes:

- Triple \(\approx 1/48\) is unchanged in all cases (weight is \(O(\varepsilon)\)).
- Overlap-fixed double \(=-0.003481\approx -1/288\); projected double residual \(\approx 0\).
- Literal global doubles the soft-soft ll residue (\(-32\)) and also shifts pl by \(-48\) relative to field-only, so \(a_1\to a_1-64\) and B double \(\approx+0.017352\neq-1/288\). Projected double cancelation is **spoiled**.

## Verdict on \(-16\,/\,64/3\)

| Reading | Result |
| --- | --- |
| Soft-soft residues of the literal vertex/coupling operator | **Confirmed:** ll \(=-16\), \(W/U=+64/3\) |
| Sector-restricted overlap assembly (\(\Delta a_1=-16\), \(\Delta\mathrm{pl}=0\), \(\Delta c_0=64/3\)) | **Confirmed:** unique among the three weight choices that restores projected double cancelation |
| Literal *global* Lagrangian doubling (\(w_f=w_c=1\)) | **Refuted** as the Newton ledger story: shifts pl, B double \(\neq-1/288\) |

What remains open is scheme/assembly: why the second operator copy is present only on soft-soft overlap and not globally — not whether \(-16\) and \(64/3\) are the soft-soft numbers of the literal vertex operator (they are).

## FFff(A) / CcFf / QFf consistency

Same identity \(\varepsilon=\lambda_F\) applies at every Yukawa vertex:

- **FFff(A):** field-weight already uses \(-\lambda_F\sum T_{ii}\) on non-root vertices; a literal global coupling copy would add another \(-\varepsilon\sum T\) and shift A poles — the present ledger does not.
- **CcFf:** one Yukawa pair \((F,f)\) plus cubics; same \(\varepsilon\) per Yukawa vertex.
- **QFf:** quartic tadpole + Yukawa sunset; Yukawa weight \(\varepsilon\) sits in the sunset.

Whatever convention is chosen for B (field-only / field+coup global / overlap-only) must be applied uniformly to A, CcFf, and QFf. The ledger’s consistency with projected double cancelation selects the overlap-restricted (or equivalently field-only plus soft-soft-only addback) convention for all of them.
