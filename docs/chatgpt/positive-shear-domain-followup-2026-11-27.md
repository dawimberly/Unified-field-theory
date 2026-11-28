# Positive shear and the clock domain: construction attempt

November 27, 2026. Follow-up to [clock-domain-audit-2026-11-27.md](clock-domain-audit-2026-11-27.md).

**Result:** I did not construct a dense domain for the full nonlinear ordered gravitational clock. I found a concrete obstruction in a free null-mode benchmark: retaining the unrenormalized positive shear contribution makes the associated exponential clock diverge for every state. A subtraction repairs that benchmark, but its extension to the full clock is not justified. This is a narrower result than a counterexample to the paper.

## Source grounding and exact scope

The target is Chandrasekaran, arXiv:2611.34924v1, equations (4.22), (4.24)–(4.27), https://arxiv.org/pdf/2611.34924. The focusing integrand is positively ordered; the clock sandwiches its exponential between area factors. Its canonical fields are functions on a reduced null phase space, so they cannot automatically be treated as jointly commuting configuration variables. The canonical embedding and the area factors must be retained in any full proof.

Before new tests I read the relevant quantization/ordering passages of two peer-reviewed primary papers: Ciambelli–Freidel–Leigh, *JHEP* **12** (2024) 028, §3, https://doi.org/10.1007/JHEP12(2024)028; and Hollands–Wald, *Commun. Math. Phys.* **223** (2001) 289, introduction and construction of Wick fields, https://arxiv.org/pdf/gr-qc/0103074. They establish the relevance of ordering and composite-field prescriptions; neither provides the missing full-clock domain proof.

Everything below is my independent benchmark calculation. It is not an independent reproduction of the paper's nonlinear operators.

## 1. Construction target

For increasing clock forms Q_c, the earlier exact variational inequality was

\[
1-\langle\psi,L(L+T_c)^{-1}\psi\rangle
\leq \|\psi-\eta\|^2+Q_c[\eta]/L.
\]

Thus we need trial vectors approaching the reference state while their actual clock forms remain bounded uniformly in the cutoff. A cutoff-dependent normalizable state is insufficient if it has no nonzero norm limit in the chosen representation.

There are two separate limits: a focusing/exponential cap c and a radial field-mode cutoff N. The preceding audit concerned c. Uniformity in N is an additional continuum requirement, not something supplied by a fixed-resolution argument.

## 2. Why a commuting configuration model is insufficient

For a free real null profile, the radiative symplectic form has the structure

\[
\Omega\propto\int\delta(\partial_\theta z)\wedge\delta z\,d\theta.
\]

Writing a Fourier mode as a cosine coefficient plus a sine coefficient pairs those coefficients symplectically. They therefore become conjugate quadratures. Multiplication by smooth compactly supported functions of *all* profile values would remove this noncommutativity and would solve a different problem. I do not use that shortcut to claim a domain for the original clock.

## 3. A null-mode benchmark retaining both quadratures

Define a real periodic, single-polarization free profile on a dimensionless circle. Let s>0, gamma>0 and

\[
z_N(\theta)=\sum_{k=1}^N\frac{Q_k\cos k\theta+P_k\sin k\theta}{\sqrt{s k}},
\qquad[Q_k,P_l]=i\delta_{kl}.
\]

The constants absorb symplectic normalization. This is a deliberately specified benchmark: periodic domain, one polarization, uniform area, constant smearing, quadratic focusing. It replaces the full partial-null-cone focusing integral and omits nonlinear and area-density operators.

Fourier orthogonality gives the exact positive quadratic operator

\[
X_N=\frac{\gamma}{\pi}\int_0^{2\pi}(\partial_\theta z_N)^2d\theta
=\frac{\gamma}{s}\sum_{k=1}^Nk(Q_k^2+P_k^2)
=\frac{\gamma}{s}\sum_{k=1}^Nk(2n_k+1).
\]

The Gaussian vacuum has profile cosine variance 1/(2sk), giving the standard 1/k free-mode scaling. Unlike a coordinate-only model, the conjugate sine coefficient is retained.

Since n_k>=0,

\[
X_N\geq c_N I,\qquad c_N=\frac{\gamma N(N+1)}{2s}.
\]

For the **benchmark** clock T_N=exp(X_N), the spectral theorem gives

\[
T_N\geq e^{c_N}I,\qquad
\|L(L+T_N)^{-1}\|\leq\frac{L}{L+e^{c_N}}\longrightarrow0
\]

for every fixed L as N grows. This is an operator-norm statement, not merely a divergent vacuum mean. No normalized trial state, squeezed state or finite-particle superposition can have a uniformly bounded T_N form. The increasing limit form has only the zero vector in its finite domain.

If a separate spectral cap c is imposed on X_N, then T_(N,c)=exp(min(X_N,c)) is bounded at fixed c. Removing that cap first restores T_N and the obstruction. At fixed c, taking N to infinity instead gives exp(c)I; subsequently c→infinity again loses every finite reading. A finite cap is not a continuum repair here.

For s=100 and gamma=1:

| N | c_N | Upper bound on resolvent norm |
| ---: | ---: | ---: |
| 8 | 0.36 | 0.8746 |
| 16 | 1.36 | 0.7196 |
| 32 | 5.28 | 0.04846 |
| 64 | 20.8 | 9.26e-9 |
| 128 | 82.56 | 1.40e-35 |

### Why squeezing is not enough

For a centered pure squeezed oscillator,

\[
\langle Q_k^2+P_k^2\rangle=\cosh(2r_k)\geq1.
\]

The lower bound is state-independent by the oscillator uncertainty relation. Shrinking one profile quadrature enlarges its conjugate. State selection cannot subtract this benchmark's vacuum term.

## 4. A second attempted construction: damping high modes

As a separate **commuting** diagnostic, take independent q_k of variance 1/(2sk), and S_N=sum k^2 q_k^2. Try

\[
\eta_N=\frac{e^{-\beta S_N}\Omega}{\|e^{-\beta S_N}\Omega\|}.
\]

For beta>0 the exact overlap is

\[
\langle\Omega,\eta_N\rangle
=\prod_{k=1}^N\frac{(1+2\beta k/s)^{1/4}}{(1+\beta k/s)^{1/2}}.
\]

Its logarithm is -N log(N)/4+O(N), so it tends to zero. Moreover, eta_N is not Cauchy: its overlap with eta_(2N) is the same product over the newly damped modes and tends to zero. The norm distance approaches sqrt(2).

For the exponential moment of exp(c S_N) to exist, every retained mode requires

\[
c<2\beta+s/k.
\]

At fixed positive c, beta therefore cannot tend to zero while all high modes are retained. Making a fixed-cutoff damped state look arbitrarily close to the vacuum is incompatible with this construction's high-mode requirement.

At s=100 and beta=0.5, vacuum overlaps for N=32,128,512,2048 were approximately 0.9443, 0.1139, 4.05e-21, and 3.19e-244. This rejects that damping family as approximants in the original infinite tensor-product sector. It is not a theorem about all possible gravitational states.

## 5. A repair that works in the benchmark

Subtracting the vacuum constant yields

\[
X_N^{\rm ren}=X_N-c_N I
=\frac{2\gamma}{s}\sum_{k=1}^Nk n_k.
\]

On the bosonic Fock space over countably many modes, define the positive diagonal operator Xren by this occupation sum. Finite linear combinations of finite-occupation basis vectors form a dense domain for exp(Xren/2), and in fact a core by spectral truncation. For any such vector, its exponential-clock form becomes independent of N once its occupied modes have been included. The vacuum has Xren=0. Thus the earlier variational recovery proof now applies on this benchmark's dense set.

This explicitly demonstrates a successful domain construction after the specified subtraction. It changes the operator, not just the state.

### Why this is not yet a gravitational repair

The full-period free oscillator energy remains nonnegative after subtraction. A local square generally does not: already Q^2-<0|Q^2|0>=Q^2-1/2 has negative expectation in a sufficiently Q-squeezed state. A finite projected oscillator check exhibited a negative eigenvalue near -0.355.

Consequently, replacing every local positive square in the original expression by a Wick square need not preserve the nonnegative focusing property used in the clock bounds. Nor do perturbative Wick-polynomial results automatically construct exponentials of the resulting nonlinear operators. Normal-ordering an exponential is another prescription and cannot be interchanged with exponentiating a normal-ordered quadratic without analysis.

I have not proved that the full geometric operators have the benchmark's state-independent divergent lower bound. The original variable weight, partial integration, angular dependence, constrained canonical embedding, and noncommuting area factors can all affect that issue. In particular, an integrated quadratic lower bound cannot simply be pushed through exp for arbitrary noncommuting operators: the exponential is not generally operator monotone. Only the scalar spectral lower bound used in the explicit benchmark justified that step.

## 6. Branch decision and next mathematical obligation

| Attempt | Bounded conclusion |
| --- | --- |
| Smooth compact support in all profile values | Does not address the noncommuting null data. |
| Bare positive periodic shear clock | Closed: exact state-independent obstruction as N→infinity. |
| Gaussian high-mode damping above | Closed: the family does not converge in the reference sector. |
| Vacuum-subtracted periodic shear clock | Successful dense domain for this free benchmark. |
| Full nonlinear ordered gravitational clock | Still open; no uniform bound or actual counterexample established. |

The next full-theory obligation is precise: specify the field-mode regulator and composite-operator prescription for (4.22), then prove that the associated forms in (4.27) admit a common dense domain with bounds surviving both cutoffs. This must retain the reduced commutators and area factors. A valid construction might differ from the free benchmark through these structures, but that mechanism must be demonstrated.

For the emergence program, this is a useful rejection criterion: a microscopic candidate must control short-wavelength vacuum fluctuations in its clock and focusing observables, while retaining the required large-scale constraints. This audit has not supplied a microscopic model or selected one.

## Reproducibility and limits

Only local CPU was used. The script below passed 38 assertions: numerical Fourier integration against the analytic positive quadratic; an oscillator matrix identity with buffer levels before projection; squeezing-energy bounds; finite-occupation stability after subtraction; and a negative local Wick-square example. Analytic statements, not the assertion count, establish the infinite-mode benchmark results. No full nonlinear gravitational simulation was run.

```python
import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln
# Periodic one-polarization free null field; dimensionless theta, [Q_k,P_l]=i delta_kl.
# z(theta)=sum (Q_k cos ktheta + P_k sin ktheta)/sqrt(s*k).
# X=gamma/pi int (dz/dtheta)^2 dtheta =gamma/s sum k(Q_k^2+P_k^2).
s=100.; gamma=1.; assertions=0
rng=np.random.default_rng(20261127)
for N in [1,2,3,8,16]:
 for _ in range(3):
  Q=rng.normal(size=N); P=rng.normal(size=N); k=np.arange(1,N+1)
  def dz(t): return np.sum(np.sqrt(k/s)*(-Q*np.sin(k*t)+P*np.cos(k*t)))
  numeric=gamma/np.pi*quad(lambda t:dz(t)**2,0,2*np.pi,epsabs=1e-10,limit=300)[0]
  analytic=gamma/s*np.sum(k*(Q**2+P**2))
  assert abs(numeric-analytic)<1e-9; assertions+=1
# Independent oscillator matrix check with buffer levels before projecting squares.
D=8; a=np.diag(np.sqrt(np.arange(1,D)),1); Q=(a+a.T)/np.sqrt(2); P=(a-a.T)/(1j*np.sqrt(2))
assert np.max(np.abs((Q@Q+P@P)[:D-1,:D-1]-np.diag(2*np.arange(D-1)+1)))<1e-12; assertions+=1
print('N, bare minimum X_N, bare vacuum clock exp(X_N), uniform R_N(10) upper bound')
for N in [1,4,8,16,32,64,128]:
 c=gamma*N*(N+1)/(2*s)
 print(N,c,np.exp(c),10/(10+np.exp(c)))
# Squeezing cannot reduce oscillator uncertainty energy: <Q²+P²>=cosh(2r)>=1.
for r in np.linspace(-4,4,17):
 assert .5*np.exp(2*r)+.5*np.exp(-2*r)>=1-1e-14; assertions+=1
# Wick-subtracted total periodic generator: Xren=2gamma/s sum k n_k.
# A fixed finite-occupation vector has cutoff-independent exponential-clock energy.
occupations={1:2,3:1,7:1}
for N in [7,16,64,128]:
 e=2*gamma/s*sum(k*n for k,n in occupations.items() if k<=N)
 assert abs(e-.24)<1e-14; assertions+=1
# Commuting coordinate-only benchmark: Gaussian damping requires a non-small UV squeeze.
# eta_N ∝ exp[-beta sum k²q_k²] Omega_N; log overlap computed exactly.
beta=.5
for N in [8,32,128,512,2048]:
 k=np.arange(1,N+1)
 log_overlap=np.sum(.25*np.log1p(2*beta*k/s)-.5*np.log1p(beta*k/s))
 print('damped-state overlap',N,np.exp(log_overlap))
# Normal ordering of a local positive square loses positivity already for one oscillator.
local_ren=(Q@Q-.5*np.eye(D))[:D-1,:D-1]
assert np.linalg.eigvalsh(local_ren)[0]<0; assertions+=1
print('local Wick-square smallest finite projected eigenvalue',np.linalg.eigvalsh(local_ren)[0])
print('PASS',assertions,'assertions; all computations are benchmarks, not the full nonlinear clock.')

```