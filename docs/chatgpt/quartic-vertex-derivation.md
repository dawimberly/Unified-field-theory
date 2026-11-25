# Independent quartic potential: tree-level expansion check

**Conversation date range** 2026-10 (exact conversation dates UNVERIFIED)

**DERIVED** This checkpoint records independent quartic potential: tree-level expansion check. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** This is a classical algebra calculation in D=4 with mostly-minus signature. It uses a simplified tree-level clock/background coefficient V0, before specifying quantum normal-ordering and counterterms. It is not a complete quantum implementation of the constrained composite theory.

**DERIVED** Let K_mu_nu=sum_a partial_mu phi_a partial_nu phi_a, A=eta_inverse K, U2=m0^2 phi^2/2, U4=lambda(phi^2)^2/(4N). Scale both phi and its derivatives by a real bookkeeping amplitude t.

**DERIVED** On the branch continuous from t=0, the determinant action density is

**DERIVED** $L(t)=V0 sqrt(det[I+t^2 A/V0]) / [1+t^2 U2/V0+t^4 U4/V0].$

**DERIVED** Define

- **DERIVED** a1=tr(A)/2;
- **DERIVED** a2=[tr(A)]^2/8-tr(A^2)/4;
- **DERIVED** a3=[tr(A)]^3/48-tr(A)tr(A^2)/8+tr(A^3)/6.

**DERIVED** Multiplying the determinant and reciprocal-potential series gives

**DERIVED** $L(t)=V0+t^2 L2+t^4 L4+t^6 L6+O(t^8),$

**DERIVED** $L2=a1-U2,$

**DERIVED** $L4=(a2-a1 U2+U2^2)/V0-U4,$

**DERIVED** $L6=(a3-a2 U2+a1 U2^2-U2^3)/V0^2 + (2U2-a1)U4/V0.$

**DERIVED** Thus the independent quartic potential contributes the familiar -U4 at four-field order AND new six-field terms

**DERIVED** Delta L6 = [2 U2 U4 - (1/2)tr(A) U4]/V0.

**DERIVED** These include both a potential term and a derivative term. Keeping -U4 but discarding Delta L6 does not reproduce the exact composite action through sixth order. Counterterm choices and the field expansion about the selected quantum background must be worked out before assigning the full loop-diagram vertices.

## Independent checks

**DERIVED** Twelve random six-component fields and derivative matrices were generated with a fixed seed. The exact determinant was evaluated by its 24-term permutation formula at 60-digit precision, independently of the trace formula used for the truncated series. The amplitudes are t=0.08,0.06,0.04,0.03,0.02. Double precision was insufficient for reliable asymptotic subtraction in some cases, which motivated this precision choice.

**DERIVED** The correct remainder slopes range from 7.9683 to 8.0217. Omitting Delta L6 gives slopes from 5.9996 to 6.0019. A separate derivative with respect to lambda agrees with the new sixth-order coefficient. All 96 domain, convergence, omitted-term and derivative checks pass. These checks validate this expansion only. They do not establish a massless graviton, conservation at quantum loop order, or originality of the algebra.

**DERIVED** Run `python quartic_composite_vertex_audit.py` to reproduce results.json. The original determinant-action framework and its expansion are in the primary papers listed in RESEARCH_BRIEF.md; the quartic extension here is elementary bookkeeping based on that action.
