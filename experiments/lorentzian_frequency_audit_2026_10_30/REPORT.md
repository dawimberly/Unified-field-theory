# Conditional Lorentzian frequency and causal-response audit

Date: 2026-10-30 US Central. Baseline: 945106beb3436dca4eb039926eb9686b635b937e.

This batch applies an explicit frequency continuation to the complete saved first-sheet coefficient pole inventories. It selects neither a microscopic metric integration cycle nor a constrained quantum theory. The point is to distinguish three claims that must not be merged: temporal support after a source, absence of exponential growth, and a positive/physical quantum response.

For the exact formally matched kernels at k=0, the numbers of additional upper-frequency-half-plane poles are 3,3,0,5,3 for L=3,5,10,30,100, summed across scalar and spin-2 coefficients without counting tensor polarization multiplicity. These are images of existing Euclidean poles, not newly discovered modes. Raw kernels have no such additional upper poles; their known negative scalar continuum regions still prevent interpreting this absence as physical health.

## Literature-access limitation

The primary-literature retrieval service returned empty error results for both a targeted search and direct author-manuscript requests. No new paper was read and no detailed result from those inaccessible articles is claimed. sources.json records the attempted URLs strictly as unverified retrieval targets. The derivations below are elementary frequency-map and residue calculations with their assumptions explicit. A peer-reviewed treatment of runaway-removal prescriptions and effective-theory order reduction remains pending. This temporary access problem did not block the mathematical batch or justify selecting a physical contour.

## 1. Assumptions and branch map

Use e^(-i omega t) for inverse time Fourier transforms, real spatial momentum magnitude k, and the formal substitution

q=k^2-(omega+i0)^2.

The Euclidean coefficient inverses are F0(q)=1/A(q) and F2(q)=1/K2(q), with their saved principal parameter logarithms and cut q<=-4. This is a specified analytic continuation of those functions. It is not assumed automatically equal to a retarded commutator in an unspecified quantum state or constrained measure.

For omega=u+iv, v>0, Im q=-2uv. If Im q=0 inside that upper half-plane, then u=0 and q=k^2+v^2>0. Thus the upper frequency half-plane never maps to the Euclidean negative-real cut. Every constituent/PV threshold maps to the real frequency boundary omega=+/-sqrt(k^2+4z_i). The i0 selects the lower q lip for positive real omega and the upper lip for negative omega. The regulator cuts do not create a new upper-frequency branch cut in this substitution.

A simple Euclidean pole at q_j maps to

omega_j=+/-sqrt(k^2-q_j),
r_omega=-r_q/(2 omega_j),

away from omega_j=0. Both square-root solutions are retained; picking one by its sign would discard part of the response. The infinitesimal i0 shifts the frequency poles by -i0 and does not remove a finite positive imaginary part.

- q_j=-s<0 gives real boundary frequencies +/-sqrt(k^2+s). These locations alone do not cause exponential growth.
- q_j=a>0 gives one positive imaginary frequency when k^2<a, with growth rate sqrt(a-k^2). For k^2>a the two frequencies are real. At k^2=a the frequency zero is double and needs a separate limiting treatment.
- A nonreal q_j gives exactly one upper and one lower frequency pole for every finite k. A conjugate q pair therefore gives two upper and two lower poles. Reality pairs the upper poles as omega and -conjugate(omega).

For a complex q_j=x+iy, the positive growth rate is

gamma(k)=sqrt([sqrt((k^2-x)^2+y^2)-(k^2-x)]/2)>0.

At large k it approaches |y|/(2k), while the real oscillation frequency is large. A small imaginary frequency component of a high-q pole is not a controlled low-energy prediction.

## 2. What upper poles do and do not imply

For a tempered retarded distribution with support t>=0, multiplying by exp(-v t), v>0, gives its Fourier-Laplace transform analytic in Im omega>0, with the appropriate polynomial bounds. Such a response cannot contain a genuine pole in that open half-plane. This observation does not require a positive spectral density.

Therefore, if a source has nonzero overlap with one of the mapped upper poles, the exact specified frequency function cannot at the same time be the transform of a tempered, non-exponentially-growing causal response. A local source contact changes only distributions at t=0 and cannot remove an isolated upper pole. An analytic remainder cannot cancel its residue unless the complete function actually has no pole there.

Causal support by itself is weaker. A causal inverse-Laplace response can grow exponentially, and its transform then need only be analytic above a positive growth abscissa rather than throughout Im omega>0. Calling every upper pole 'acausality' would confuse causal growth with stable response.

Conversely, choosing a frequency contour that removes future runaway growth can change temporal support. Such a choice is extra boundary/contour information, not proof of a microscopic physical prescription. No Lee-Wick, future-boundary, retarded or other metric path-integral contour is adopted here.

## 3. Explicit pole-sector contour calculation

Isolate one already saved coefficient pole and write

f(omega)=r_q/(k^2-omega^2-q_j)=r_q/(Omega^2-omega^2),
Omega^2=k^2-q_j.

An inverse-Laplace line above both frequency poles gives, by closing below for t>0 and above for t<0,

g_causal(t)=theta(t) r_q sin(Omega t)/Omega.

It obeys g(0+)=0 and g'(0+)=r_q. When Im Omega is nonzero, its positive-time behavior generically contains an exponentially growing term. Conjugate q sectors can be combined to obtain real test responses; an isolated complex component is not itself a real physical pulse.

For the same pole sector with no real-axis frequency poles, an ordinary real-frequency contour instead gives

g_real(t)=i r_q exp[-i Omega_minus |t|]/(2 Omega_minus),
Im Omega_minus<0.

It decays as |t| increases but is generally nonzero at negative t. For a positive Euclidean pole q_j=a and k^2<a, put gamma=sqrt(a-k^2). The formulas reduce to

g_causal=theta(t) r_q sinh(gamma t)/gamma,
g_real=-r_q exp(-gamma |t|)/(2gamma).

These are different inverses of the same differential pole factor with different boundary conditions. A damping-looking real-axis inverse is not a causal stability proof. The examples isolate pole terms only; no full-kernel time evolution, nonlinear stability, spatial propagation/light-cone theorem or quantum contour construction was executed.

## 4. Apply the complete inventory, preserving source distinctions

At k=0 the extra upper-frequency counts are:

| L | Raw scalar/spin-2 | Matched scalar | Matched spin-2 | Matched total |
|---|---:|---:|---:|---:|
| 3 | 0 / 0 | 2 | 1 | 3 |
| 5 | 0 / 0 | 2 | 1 | 3 |
| 10 | 0 / 0 | 0 | 0 | 0 |
| 30 | 0 / 0 | 3 | 2 | 5 |
| 100 | 0 / 0 | 1 | 2 | 3 |

The massless origin is omitted. Its k=0 double-frequency zero and usual retarded secular limit are separate from the additional exponential poles. Matched L=10 is the unchanged reference. The raw subthreshold scalar poles map to real frequency boundary pairs, not extra upper poles.

For example, at k=0 the L=30 collar pair has gamma approximately 0.303666, whereas its matched spin-2 pair has gamma approximately 6.60584, in units with constituent mass one. These are numerical continuation diagnostics, not rates in seconds or predictions for a physical vacuum. All extra unstable mapped poles originate at Euclidean |q|>4, outside the low-energy Taylor radius. Matching finitely many Einstein/curvature coefficients does not require a UV completion to reproduce this full high-q analytic function.

The previously fixed constituent scalar source multiplies a pole residue by (q_j-2)^2/16. The known additional poles have nonzero overlap factors in that channel; the kinematic zero at q=2 removes none of them. Spin-2 source overlap must be treated separately. A generic conserved transverse-traceless source can overlap a spin-2 pole. Different tensor channels, massless gauge constraints and source-dependent cancellations cannot be replaced by adding scalar signs in isolation. No additional particle count or S-matrix conclusion follows from the table.

An exact analytic response agreeing with the continued function on an open interval of interior imaginary frequencies is fixed by analytic continuation where connected; deleting its upper poles is not an innocuous unchanged analytic response. Agreement of only a finite low-energy expansion is much weaker and allows different UV functions. Thus this conditional obstruction neither rules out every UV completion nor transfers to a constrained theory whose source response has not been matched.

## 5. New execution evidence

Run `python experiments/lorentzian_frequency_audit_2026_10_30/check.py`.

The script reads saved pole inventories, residues and spin-2 certificate centers. It does not search for or recertify Euclidean roots. It executes:

- 88 frequency-map entries at k=0,1,10,100, preserving both frequency roots and both members of every conjugate pair.
- 56 symmetric near-frequency-pole limits checking the derivative-chain residue sign and normalization; maximum normalized discrepancy 1.80e-10.
- 12 independent oscillatory real-frequency quadratures against the closed two-sided inverse, including negative times; maximum normalized discrepancy 6.37e-13.
- Causal jump checks and saved positive/negative-time pole-sector examples.

The frequency existence/sign classification follows algebraically from the certified q locations, not from these floating-point samples. Displayed numerical frequencies and residues are not interval-certified digits. The numerical checks normalize by max(1,|expected|). No previous successful Ward, contour count, source dispersion or broad pole scan was repeated; no assertion failed or tolerance was relaxed.

The present conclusion is conditional and limited: the four non-reference matched prescriptions have upper frequency poles under the explicit continuation, so a source-overlapping response cannot simultaneously have this function and be tempered causal. A causal growing pole-sector inverse and a decaying two-sided inverse illustrate the distinction. Whether a microscopic theory selects either, modifies the response, or uses a different state/source definition remains open. Absence of such poles in the raw prescriptions is not a unitarity or health proof.

Next finish the primary-literature comparison of runaway-removal prescriptions and effective-theory order reduction when retrieval returns. Then prepare a concrete research decision memo identifying which unresolved microscopic source/measure or regulator assumption must be specified to move beyond this fully audited leading-N auxiliary calculation. Do not spend the next batch repeating already complete pole inventories or present a contour choice as a routine calculation.
