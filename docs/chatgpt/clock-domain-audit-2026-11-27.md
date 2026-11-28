# Clock-domain audit — November 27, 2026

**Verdict:** a general shortcut to the clock-domain claim is obstructed. A gravity-specific proof remains open; this audit does not disprove the proposed theory.

## 1. Exact point under examination

Source: V. Chandrasekaran, *An algebraic quantization of gravity from the null cone*, arXiv:2611.34924v1, https://arxiv.org/pdf/2611.34924.

Section 5.1, footnote 18, labels the passage from a classical focusing-tail estimate to the exact noncommuting matrix element heuristic. Appendix C.2, (C.13)–(C.17), explicitly retains the projection onto the closure of the limiting form domain; its coherent-state density argument is conditional on the missing estimate. The classical reasoning is at fixed resolution. Neither a formal expansion nor that classical estimate alone establishes the required exact quantum limit.

For mathematical grounding, I read the relevant monotone-form results and proofs in J. Behrndt, S. Hassi, H. de Snoo and R. Wietsma, *Proceedings of the Royal Society of Edinburgh A* **140** (2010), 927–951, especially Theorems 1.3 and 4.2. DOI: 10.1017/S030821050900078X. Author manuscript: https://www.numerik.math.tu-graz.ac.at/~behrndt/BehrndtHassiSnooWietsma09.pdf. Those theorems permit a nondensely defined limiting form and a self-adjoint relation, rather than necessarily an operator on the original Hilbert space.

## 2. Independent derivation: what must actually be proved

Work on a fixed-area reference sector H (the paper's K_A), not automatically the entire extended representation. Let Q_c be increasing nonnegative closed forms, with positive operators T_c. Assume their monotone limit Q is closed, as supplied by the monotone-form theorem. Put M=closure(dom Q), with orthogonal projection P_M; T denotes the operator associated with Q on M.

For L>0 define R_c(L)=L(L+T_c)^(-1). The strong limit is

    R_infinity(L)=L(L+T)^(-1) P_M.

The spectral theorem on M gives

    lim_(L→infinity) R_infinity(L) psi = P_M psi.

Consequently, for a normalized psi,

    lim_(L→infinity) [1-<psi,R_infinity(L)psi>]
       = ||(I-P_M)psi||^2.

This is an exact diagnostic: the missing estimate is equivalent to proving psi belongs to M. It cannot be obtained by silently treating the limiting clock as an everywhere-supported self-adjoint operator; that would assume the domain result.

An infinite clock expectation is not itself a failure. A vector can lie in M while outside dom(Q), and its resolvent approximation can still converge in norm.

### A sufficient replacement bound

For any normalized psi and trial vector eta in dom(Q_c), completing the square gives

    d_c(L) := 1-<psi,R_c(L)psi>
       = inf_eta { ||psi-eta||^2 + Q_c[eta]/L }.

The minimizer is R_c(L)psi. Therefore, if for every epsilon>0 there is one eta_epsilon such that

    ||psi-eta_epsilon|| < epsilon,
    sup_c Q_c[eta_epsilon] <= C_epsilon < infinity,

then

    sup_c d_c(L) <= epsilon^2 + C_epsilon/L,

and the desired iterated limit is zero. For a dense collection of psi, this establishes M=H. A dense set of uniformly bounded-form trial vectors suffices; it need not be shown to be a form core in the stronger graph-norm sense.

This replaces an uncontrolled classical-to-quantum comparison with an exact variational inequality. The missing gravitational work is to construct these approximants for the actual ordered clock, including its area weights. Merely cutting off each noncommuting factor independently does not establish a common dense domain.

### Equivalent spectral target

Let mu_c,psi be the spectral probability measure of T_c in psi. Then

    d_c(L) = integral u/(L+u) dmu_c,psi(u).

For any K>0,

    (1/2) mu_c,psi([L,infinity)) <= d_c(L)
       <= K/L + mu_c,psi((K,infinity)).

Thus uniform spectral tightness, sup_c mu_c,psi((K,infinity))→0 as K→infinity, is another sufficient and, here, equivalent condition for uniform resolvent recovery. This is a statement about quantum spectral measures; a positive classical Wigner density does not supply it without an additional comparison theorem.

## 3. Independent counterexamples and controls

These examples test logical sufficiency. They are not the gravitational clock and do not reproduce all its structural assumptions.

### Positive monotone clocks can lose a subspace

On C^2 let P=diag(1,0), T_n=I+(n-1)P, and psi=(1,1)/sqrt(2), n>=1. Every T_n is positive and bounded, and the family increases. Explicitly,

    d_n(L) = n/[2(L+n)] + 1/[2(L+1)].

Removing n first yields

    lim_n d_n(L) = 1/2 + 1/[2(L+1)],
    lim_L lim_n d_n(L) = 1/2.

The limit form domain is ker(P), not dense. Reversing the limits gives zero. Along n=sqrt(L), d_n(L) also tends to zero. This shows why control along a growing finite cutoff does not by itself control removal of the cutoff first. With n=exp(c), these clocks even obey T_0<=T_c<=exp(c)T_0.

Taking T_n=nI instead makes M={0}. Positivity, finite-cutoff boundedness and monotonicity alone therefore do not even guarantee a nonzero finite-reading sector.

### Noncommutativity does not remove the obstruction

Let B=[[1,0.3],[0.3,1]] and T_n=B[I+(n-1)P]B. These are again positive monotone clocks, bounded between B^2 and nB^2. Distinct T_n do not commute. The limit domain is ker(PB), a proper closed line. If w is its unit vector, then

    R_infinity(L) = L/[L+<w,B^2 w>] |w><w|.

Thus recovery is again a projection, not the identity. This rules out the hope that merely adding noncommutativity repairs the general implication. It does not rule out a special gravitational estimate.

### Every perturbative coefficient can miss the loss

In the first example choose psi_h=(sqrt(p_h),sqrt(1-p_h)), p_h=exp(-1/h^2), h>0. Then

    lim_n d_n(L) = p_h + (1-p_h)/(L+1),
    lim_L lim_n d_n(L) = p_h > 0.

The formal expansion about h=0 contains only 1/(L+1): every coefficient of the extra term is zero. Hence coefficientwise O(1/L) control cannot exclude a finite-h domain defect without a remainder estimate. This is an abstract beyond-all-orders example, not a proposed size or mechanism for an effect in gravity.

### Successful control

On l^2(N), let T_n e_k=min(k,n)e_k and psi_k=sqrt(3)2^(-k), k>=1. The limiting operator has a dense domain and

    sup_n <psi,T_n psi> <= 3 sum_k k/4^k = 4/3,
    sup_n d_n(L) <= 4/(3L).

The replacement bound succeeds in this case. Conversely, choosing normalized psi_k proportional to 1/k gives an infinite limiting mean but still norm recovery by dominated convergence in the spectral sum. A divergent mean alone is not grounds for rejecting the clock.

## 4. Dependency map

| Claim or calculation | Consequence of the unresolved domain estimate |
| --- | --- |
| Gaussian preparation in a finite-reading diamond, §5.1 | Not established by the audited estimate alone. |
| Dense clock form domain in K_A, C.17 | Conditional on exact recovery for nearby coherent states. |
| Finite-dimensional diamond space, §5.2 | Its nontrivial physical realization remains conditional; this audit does not refute the separate dimension argument. |
| Diamond entropy/topology conclusions, §§5.3–5.4 | Inherit the sector-existence dependency and have further assumptions. |
| Axiomatic realization for diamonds, E.3 | Inherits that dependency. |
| Classical constraints; formal semiclassical calculations; pair-overlap example | Not invalidated by this domain objection. |

This map identifies proof dependencies, not independent validation of the other arguments. Failure of recovery for one Gaussian would not establish that every finite-reading state is absent. Even M=H would not by itself prove finite dimension, representation independence, or the continuum construction.

## 5. Verification and decision

The accompanying Python check used only local CPU, NumPy and a fixed random seed. It passed 235 numerical assertions: exact diagonal resolvents; a noncommuting limiting resolvent; 100 variational identities and trial-vector inequalities; and the successful uniform-energy control. Numerical checks corroborate the analytic identities above; finite matrices cannot decide the gravitational infinite-dimensional domain question.

**Closed route:** positivity plus monotone convergence plus a finite-cutoff diagonal estimate plus coefficientwise perturbative suppression is not a general proof of full quantum-clock recovery. Explicit counterexamples obstruct that inference.

**Open route:** construct a dense set of trial vectors with uniformly controlled actual clock forms, or prove quantum spectral tightness. Neither has been established here. No counterexample to the full gravitational construction has been produced.

**Implication for emergence:** use this as a domain and cutoff-consistency requirement. It supplies no derivation of gravitational constraints from microscopic dynamics, and does not require discarding our emergence objective.

## Reproducible numerical check

Run with Python 3 and NumPy. This is a logical stress test, not a numerical discretization of gravity.

```python
"""Analytic identities are in the accompanying report; these are numerical checks."""
import numpy as np
from numpy.linalg import solve, norm, eigh
I=np.eye(2); P=np.diag([1.,0.]); psi=np.array([1.,1.])/np.sqrt(2)
def R(T,L): return L*solve(L*I+T,I)
def loss(T,L,v): return float(v@(I-R(T,L))@v)
checks=0
# Monotone positive clocks, uniformly bounded at each finite cutoff.
for n in [1.,10.,100.,1e4,1e8]:
 T=I+(n-1)*P
 for L in [2.,10.,1000.]:
  exact=.5*n/(L+n)+.5/(L+1)
  assert abs(loss(T,L,psi)-exact)<1e-12; checks+=1
# Diagonal cutoff n=sqrt(L) hides the escaped half of the state.
for L in [1e2,1e4,1e6,1e8]:
 print('diagonal',L,loss(I+(np.sqrt(L)-1)*P,L,psi),'removed-cutoff lower bound',.5)
# A congruence supplies genuinely noncommuting clocks.
B=np.array([[1.,.3],[.3,1.]]); T0=B@B
v=B@np.array([1.,0.]); w=np.array([-v[1],v[0]]); w/=norm(w)
PM=np.outer(w,w); QonM=float(w@T0@w)
assert norm(T0@(T0+B@P@B)-(T0+B@P@B)@T0)>1e-3; checks+=1
for L in [2.,20.,200.]:
 T=T0+(1e9-1)*B@P@B
 expected=L/(L+QonM)*PM
 assert norm(R(T,L)-expected)<2e-6; checks+=1
 print('noncommuting',L,'limit-state norm',norm(expected@psi),'loss',float(psi@(I-expected)@psi))
# Variational identity and recovery inequality, for arbitrary trial vectors.
rng=np.random.default_rng(20261127)
for _ in range(100):
 A=rng.normal(size=(2,2)); T=A.T@A; L=float(rng.uniform(.1,30))
 v=rng.normal(size=2); v/=norm(v); eta=rng.normal(size=2)
 z=R(T,L)@v
 objective=norm(v-z)**2+z@T@z/L
 assert abs(objective-loss(T,L,v))<2e-12
 assert loss(T,L,v)<=norm(v-eta)**2+eta@T@eta/L+1e-12; checks+=2
# Dense-domain commuting control: geometric state, finite mean <= 4/3.
k=np.arange(1,501,dtype=float); weights=3*4.**(-k)
for n in [1,10,100,1000]:
 for L in [1.,10.,100.,1000.]:
  t=np.minimum(k,n)
  val=np.sum(weights*t/(L+t))
  assert val<=4/(3*L)+1e-14; checks+=1
# Nonperturbative missing weight: every power-series coefficient at h=0 is zero.
for h in [.5,.35,.25]:
 p=np.exp(-1/h**2)
 print('beyond-all-orders h',h,'missing norm squared',p)
print('PASS',checks,'numerical assertions; counterexamples are abstract, not gravitational.')

```