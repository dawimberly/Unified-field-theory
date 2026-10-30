# Composite source and auxiliary bubble bridge: bounded audit

Date: 2026-10-30 US Central. Baseline: c24dbcfec4ca891506b10131f7ffc5cd70876207.

## Primary-literature audit (limited source synopsis)

Batz, Erlich and Mrini, *Background-Independent Composite Gravity*, arXiv:2011.06541v1, associated DOI 10.1088/1361-6382/abefbd. Read the author PDF at https://arxiv.org/pdf/2011.06541, especially pp.3,5--7,9--16. Publisher access was robots-blocked; no final-version comparison is claimed.

The metric is composite, not independently integrated (p.3). Clock terms cancel algebraically; the clock profile selects coordinates, with a noncritical-coordinate qualification (2.2--2.6). Dimensional regularization proxies an unspecified physical regulator. Vacuum consistency is studied through adiabatic order two. A stress-tensor interaction and leading bubble chain produce the long-distance exchange (3.2--3.18). Section IV.A discusses contact terms in metric variations (4.1--4.4). These statements do not establish an unrestricted auxiliary-metric quantum measure or a finite-momentum PV spectrum.

Inspection targets for follow-up: operator renormalization in (2.25--2.27), metric-variable/source conventions in (3.11) and (4.1--4.4), and the discarded higher-derivative terms after (2.22). This audit identifies assumptions; it neither reproduces the entire paper nor alleges an error.

## New calculation: a conditional algebraic bridge

The following derivation is performed here for a separately specified quadratic stress-channel approximation. It is not attributed to an exact microscopic equivalence theorem. All quantities are finite-dimensional to make the algebra explicit; continuum regularization, operator products and gauge zero modes require separate treatment.

Use an orthonormal basis of symmetric D-dimensional tensors, n=D(D+1)/2. Let v be the trace vector (D entries equal to one, remaining entries zero) and define

P=(D-2) I/2-vv^T/2,
P^-1=2I/(D-2)-vv^T/(D-2), D>2.

The traceless eigenvalue of P is (D-2)/2; its trace eigenvalue is -1. P is invertible but indefinite. Introduce V nonzero, a symmetric invertible stress-channel covariance B, and the formal Euclidean-sign quadratic expression

S_quad(T,h)=T^T B^-1 T/2+V h^T P^-1 h/4-h^T T/2.

The sign here defines the test algebra. It is not an unannounced Wick rotation of the repository's Lorentzian theory. Completing the h square gives the reduced expression

S_red(T)=T^T[B^-1-P/(2V)]T/2,

and the linear composite O=P T/V. The h Hessian is H=V P^-1/2. Define algebraic inverse matrices where nonsingular:

C_T=[B^-1-P/(2V)]^-1,
C_O=P C_T P/V^2,
D_h=[V P^-1/2-B/4]^-1.

The exact Schur identity is

D_h=C_O+2P/V.

Therefore the auxiliary and linear composite two-point functions differ by a local source contact in this approximation. Replacing one by the other without that term is not equality of generating functionals. If this Euclidean quadratic representation were made convergent on a specified cycle, integrating h with source j would yield a factor exp(j^T P j/V) multiplying the reduced composite-source generating function; for the present indefinite H this is a formal completing-square identity, not a demonstrated real Euclidean integral. Lorentzian i factors must be derived from the chosen contour.

Define g=P/(2V). The ordered stress bubble sum is

A_chain=g+g B g+g B g B g+...=(g^-1-B)^-1=D_h/4.

The series expansion itself needs a convergence domain; the inverse is its algebraic continuation. Noncommuting matrix order is retained. Its factor of four comes from the coupling -h^T T/2. It cannot be silently dropped when comparing exchange normalizations.

Since 2P/V is independent of momentum, C_O and D_h have the same nonlocal pole parts under these assumptions. Invertibility of P and finite nonzero V are essential; source contraction can still cancel a residue. A pole of the channel matrix need not survive a specified conserved external-source amplitude. At gauge zero modes the above ordinary inverses do not exist; one must fix gauge consistently or restrict to the physical source sector before applying an inverse identity.

## Relation to the existing determinant, and the remaining gap

The earlier algebraic saddle Hessian was

H_saddle=V/(D-2)[I-vv^T/2]=V P^-1/2.

Thus the source bridge has the same algebraic quadratic normalization. In D=4, P=P^-1, and the repository stationary contact matrix is p P/2. Setting V=p would identify this constant term. Also a convention B=2 times the repository vertex integral makes its bubble term -B/4, as in D_h^-1. These are conditional identifications of constants and insertions, not a measurement or a proof that the two theories choose the same B, p, renormalization or counterterms.

The missing work is now precise: derive the regulated normal-ordered stress insertion B and the constant coupling V from one specified composite generating functional, with the same local source counterterms, before comparing to the full determinant Hessian. The paper's small-curvature/adiabatic expansion alone supplies no full finite-q equality. Its contact cancellation cannot be transplanted into our PV calculation by dropping the stationary contact matrix. The earlier Ward test already showed that dropping that matrix fails in our convention; it was not rerun here.

No full nonlinear measure equality, physical Lorentzian cycle, signature construction, curvature matching beyond the chosen orders, or constrained-CEV spectrum transfer follows. The previous finite-dimensional conditional-source formula is sharpened here to an explicit tensor normalization at the leading quadratic stress-channel level. This is a bounded matching step, not a new gravity model.

## Reproducible checks

Run `python experiments/composite_source_bridge_2026_10_30/check.py` from the repository root. It imports only the previous standard-library rational matrix helper and performs eight cases: D=3,4,5,6, each at V=11,17. B is symmetric positive definite with a diagonal plus rank-one form, deliberately not commuting with P. All inverse, Schur and normalization identities are checked with exact Fraction arithmetic. Twenty-four independent rational source quadratic forms verify the local difference. Maximum algebraic residual: exactly zero. The positive B does not remove the indefinite auxiliary integration-cycle problem.

results.json stores every case, definitions and Python version. No old root scans or Ward checks were repeated. No numerical continuum/path-integral execution is claimed.

## Next bounded batch

Build the stress bubble/contact dictionary from a regulated scalar generating functional, distinguishing the fixed-field seagull, connected stress correlator and cosmological stationary term. Compare only within the already authorized free leading-N auxiliary approximation unless a specified composite measure is available. Then give a low-energy matched response comparison with an explicit remainder bound, rather than interpreting the previously found high-q polynomial-resummed poles as low-energy predictions. Global scalar poles, interactions and 1/N remain open.
