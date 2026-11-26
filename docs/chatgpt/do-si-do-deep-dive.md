# Do-Si-Do hypothesis: finite release, incidence angle, and what the data exclude

**Research record** Stages 1–54 on hyperfield-stages, with the November 25 checkpoint and the Stage 52–54 standalone audits. Repository calculations are marked **DERIVED**, estimated parameters **FITTED**, postulates and design targets **ASSUMED**, and primary-source claims **CITED**. Unresolved source locations or mappings are marked **UNVERIFIED**.

**DERIVED** Finite release survives conditionally. Independent site release gives a visibility ceiling of roughly $0.76$–$0.90\,\mathrm{as}$ on $\tau_0$ under the stated phase mapping. A shared release with common jitter has no such visibility ceiling in the shared-phase model below. The acoustic and force-noise limits apply only to a specified mass-coupled noise extension. Under that extension the allowed peak strength depends on resonance frequency and width, and arbitrarily weak nonzero coupling survives. Neither the hydrogen residual nor the cantilever excess establishes a detection. The repo supplies no law that converts a zeptosecond timing residual into a particular collapse rate, and therefore no evidence that the two residuals share a cause.

## 1. Minimal model

**ASSUMED** Treat a lock as a prepared localized configuration and the wave as its outgoing channel. Let $q$ be the calibrated opening load, in energy units. For a simple oriented directional probe set $q=A(1+\cos\theta)/2$, where $A\ge0$ is the maximum opening load supplied by the probe and $\theta$ is its angle to the local opening axis. This is a test parametrization; it is not a universal identification of intensity, photon energy, and load. For an unoriented axis an odd response is only meaningful relative to a tagged site or another polar vector.

**ASSUMED** Use a gate $g\in\{0,1\}$ initialized closed. It opens on crossing $q=Q_+$ upward, closes on crossing $q=Q_-$ downward, and retains its state in between, with $0\le Q_-<Q_+$. The thresholds define the lock's resistance to opening and the reset load. The width $Q_+-Q_-$ is a driven hysteresis observable. Closing the gate halts further release of a remaining locked fraction; it does not recreate a released particle. A reverse capture experiment would need a specified preparation protocol and reservoir.

**ASSUMED** In a fixed geometry the fraction $L$ of an initially locked ensemble remaining locked obeys
$$
\dot L=-\frac{g(t)L}{\tau_\theta},\qquad
\tau_\theta=\tau_0+\frac{d}{2}\cos\theta,\qquad
\tau_0>|d|/2 .
$$
Here $\tau_0$ is the isotropic mean release time, and $d$ is the signed differential release coefficient for two oppositely oriented sites. The four core parameters are $(Q_+,Q_-,\tau_0,d)$. The probe amplitude and angle are inputs, not fit parameters. A calibrated experimental drive can also directly supply $q(t)$ without the illustrative angular load law.

**DERIVED** With a continuously open gate and $L(0)=1$, $L(t)=e^{-t/\tau_\theta}$ and the release-time probability density is $K_\theta(t)=e^{-t/\tau_\theta}/\tau_\theta$ for $t\ge0$. Its standard deviation is $\tau_\theta$. Finite release spread is therefore already present without a fifth parameter. An independent distribution of lifetimes, an extra dead time, or correlations across releases require additional parameters. The kernel describes time after opening; it is not automatically a correction to a measured radioactive lifetime or a coherent quantum amplitude.

**ASSUMED** For equivalent tagged sites with opposite local opening axes in H$_2$, map the release times to the relative electron phase. Set $x=\cos\beta$, where $\beta$ is the photon-to-bond angle, and $\tau_\pm=\tau_0\pm d x/2$. The phase frequency is $\omega_{\rm e}=E_e/\hbar$. **CITED** In Grundmann et al. (2020), *Science* **370**, pp. 339–341, Eq. (2), $E_e$ is the electron kinetic energy and $v_{\rm ph}=E_e/p_e$, with $p_e$ the electron momentum. The Stage 52 input $E_e=735\,\mathrm{eV}$ gives $\hbar/E_e\simeq896\,\mathrm{zs}$ (**DERIVED**). Photon energy is a separate experimental setting.

**ASSUMED fork A, independent release.** Average $e^{i\omega_{\rm e}(T_+-T_-)}$ over independent exponential site times with means $\tau_\pm$. **DERIVED** Each characteristic function is $\phi_\pm=(1-i\omega_{\rm e}\tau_\pm)^{-1}$, so the fringe coherence factor is $\phi_+\phi_-^*$. Its phase divided by $\omega_{\rm e}$ gives
$$
\Delta t=\frac{\arctan(\omega_{\rm e}\tau_+)-\arctan(\omega_{\rm e}\tau_-)}{\omega_{\rm e}}
=\frac{d\cos\beta}{1+(\omega_{\rm e}\tau_0)^2}
+O(\omega_{\rm e}^2 d^3 x^3).
$$
**DERIVED** The same coherence factor multiplies the baseline visibility by $[(1+\omega_{\rm e}^2\tau_+^2)(1+\omega_{\rm e}^2\tau_-^2)]^{-1/2}$. This is exactly $F=[1+(\omega_{\rm e}\tau_0)^2]^{-1}$ at $x=0$ and equals $F$ to leading order in $dx$. The repo's approximate observed visibility of $0.50$–$0.58$, with the unsuppressed visibility at most unity, gives $\tau_0\lesssim896\,\mathrm{zs}\sqrt{1/V_{\rm obs}-1}$, or about $0.76$–$0.90\,\mathrm{as}$. This is a conditional ceiling, with the two endpoints representing different visibility benchmarks, not a confidence interval. An intrinsic baseline below unity would tighten the ceiling. Angle-matched visibility and acceptance should be used in an exact joint fit.

**ASSUMED fork B, shared release.** A truly common random delay contributes a global phase and cancels from relative coherence. To retain the same dispersive phase law, assume the two site factors are the unit-modulus responses $\phi_\pm/|\phi_\pm|$, multiplied by a common random phase. **DERIVED** This gives the same $\Delta t$ with no release-induced visibility penalty. **UNVERIFIED** The shared phase-only response is an extra dynamical assumption. Sharing a random clock alone does not imply $F$, and using the same random exponential draw with different site rescalings still produces differential jitter. A purely shared additive delay with simultaneous responses would yield no differential signal.

**DERIVED** For either stated phase mapping, the timing fit tests $t_{\rm obs}=Rx/c+t_{\rm nd}+\Delta t$, where $R$ is the full bond length, $c$ is the speed of light, and $t_{\rm nd}$ is the standard nondipole and correlation term used in the comparison. The leading fitted odd coefficient is $d_1=dF$. A differential even term $d_2P_2(x)$, with $P_2(x)=(3x^2-1)/2$, is a separate check.

**ASSUMED** To include frequency-selective excess strength in the smallest *testable* extension, add $(\lambda_*,f_0,b)$ and use the repo's Lorentzian spectral convention
$$
\lambda_{\rm eff}(f)=\lambda_*\,\bar g\,\frac{1}{1+((f-f_0)/b)^2},
\qquad f=\omega/(2\pi),\quad Q_{\rm res}=f_0/(2b).
$$
$\lambda_*\ge0$ is the active peak spectral rate in $\mathrm{s}^{-1}$ in the repo's mass-proportional CSL yardstick per reference nucleon, $f_0>0$ is the center in Hz, $b>0$ the half-width at half maximum in Hz, and $\bar g$ is the gate exposure during the measurement. Use $\bar g=1$ for the noise benchmark or determine it from a logged, calibrated protocol. For spontaneous noise measurements this is a convention for an active channel, not evidence the laboratory apparatus opened a microscopic gate. Keep the repo's reference smearing $r_C=100\,\mathrm{nm}$ and reference mass $m_0=1\,\mathrm u$ fixed. An anisotropic force-noise coupling, amplitude-dependent equilibrium energy, and a detailed-balance violation are separate extensions.

**DERIVED** The core parameters do not fix $f_0$, $b$, or $\lambda_*$. The Fourier transform of the exponential probability kernel is a zero-centered low-pass response with squared modulus $1/[1+(2\pi f\tau_\theta)^2]$ under that transform convention. It cannot supply a finite center frequency. Equating $b=(2\pi\tau_0)^{-1}$ or $f_0=(2\pi\tau_0)^{-1}$ without a dynamics and coupling derivation would add an assumption. For instance $100\,\mathrm{zs}$ would give a low-pass half-power frequency of about $1.6\times10^{18}\,\mathrm{Hz}$, not the kHz cantilever center. The gate can exhibit a driven loop with no shift of equilibrium occupation once the drive and reset history are removed.

## 2. Constraint ledger

**DERIVED** The [stage ledger](../hyperfield/HYPERFIELD_RESULTS.md) contains Stages 1–38 and 40–53 as numbered sections; there is no heading for Stage 39. The only standalone STAGE markdown files in that directory on the reviewed branch are [52](../hyperfield/STAGE52_grundmann_504_audit.md), [53](../hyperfield/STAGE53_lpf_highfreq.md), and [54](../hyperfield/STAGE54_8kHz_4p6MHz_gap.md). The [November 25 summary](../hyperfield/SUMMARY_11-25.md) predates the latest coverage updates. The [research index](INDEX.md) and [incidence-angle notes](angle-of-incidence.md) give the earlier geometry and coupling context. Earlier hyperfield stages constrain distinct model additions: geometry and brane couplings (1–10), gravity-entanglement feasibility and classical–quantum back-action (11–29), and noise-induced bias in measurements of $G$ (30). Their limits cannot be converted into $d$, $Q_\pm$, or a spectral rate without a transfer law. In particular the Stage 28 white-floor bound $\Gamma V\le5.5\times10^{-35}\,\mathrm{m^5\,s^{-3}}$ is a CQ volume-law bound. Stage 30's verified atom case allows up to $0.057$ ppm in the repo's bias estimate, well short of the roughly $500$ ppm spread. Neither bound is applied as an independent CSL constraint.

### Timing, angular laws, gate, spread and balance

| Status and repo stage | Constraint | What is excluded, conditional on the indicated mapping |
|---|---|---|
| **FITTED**, 31/32/52 | H$_2$ at $800\,\mathrm{eV}$: matched slope $441\pm35.8_{\rm stat}\pm32_{\rm fit}\,\mathrm{zs}$, excess $+171\pm130\,\mathrm{zs}$ over the approximately $270\,\mathrm{zs}$ nondipole comparison, $1.3\sigma$. **DERIVED** repo limits $\lvert d_1\rvert<425\,\mathrm{zs}$ matched and $<489\,\mathrm{zs}$ conservative, at reported 95%; $\lvert d_2\rvert<102\,\mathrm{zs}$ from the earlier even fit. | Excludes larger additive odd slope in this observable under those systematics. The bound is on the fitted $d_1$, not directly on $d$. The independent visibility model additionally limits $\tau_0$; the shared-phase model does not. |
| **DERIVED**, 38/50/51 | He RABBITT even sensitivity $\lvert d_2\rvert\lesssim15\,\mathrm{as}$; independent He odd sensitivity floor about $28{,}700\,\mathrm{zs}$ at the reported 95% scale, without a validated conversion of the dipole–quadrupole measurement to $d_1$. Same-event H$_2$ reanalysis gives no independent likelihood. Public H$_2$ visibility contrast is about $0.074\pm0.018$ in the paired analysis. | He offers no stronger odd bound. A visibility difference is no measurement of $d_1$. Under equivalent-site relabeling, the phase shift can reverse while the fringe contrast remains unchanged. Changing the molecular-axis label rather than assuming a rotation about the beam reverses $x$: a rotation about the beam itself preserves the dot product with the beam. Detector acceptance and coincidence cuts must be tested explicitly. |
| **DERIVED**, 31 | Generic Malus-law deviation at roughly $10^{-13}\,\mathrm s$ is constrained at the $\sim1\%$ scale for the repo's tested shape; crossed X-ray polarizer leak $<2.4\times10^{-10}$ at $6.457\,\mathrm{keV}$; fixed-star energy anisotropy scale $4.6\times10^{-17}\,\mathrm{eV}$, fractional $10^{-18}$. | These concern polarization probabilities or a preferred-frame Hamiltonian, not the molecular-relative phase-delay law. The crossed-null limit only constrains deviations that survive at the null. The polarization results do not measure the duration of all releases. |
| **DERIVED**, 31 | Single microwave photons at $24.86\,\mu\mathrm{eV}$ give nondemolition registrations and linear response at sub-single-photon average population. | Rules out a universal minimum deposited energy as a prerequisite for those registrations, or a threshold above the single-quantum energy under an incident-energy identification. Mean population is not the local amplitude of an individual click; this does not calibrate $Q_\pm$ of a separate channel. |
| **DERIVED**, 31 | Population-dependent clock energy coefficient $\lvert\epsilon\rvert<2.4\times10^{-20}\,\mathrm{eV}$, corresponding to $5.8\,\mu\mathrm{Hz}$. | Constrains the projected Weinberg-like population term; an external drive acting through a gate is a different observable. |
| **DERIVED**, 31 | VIP-2 $\beta^2/2\le8.6\times10^{-31}$ at reported 90%; He direct versus shake-up extra delay $<1.5\,\mathrm{as}$ at reported 95%; detector jitter about $2.7\pm0.2\,\mathrm{ps}$. | The first applies to an exchange-distinguishable intrinsic label only with an order-one mapping; the second limits channel differences; the third is the full measured detector scale, not a calibrated bound on universal intrinsic spread. |
| **DERIVED**, 44 | $^{56}$Mn static rate fractional spread $s<0.019$, Weibull shape $\lvert\beta_{\rm W}-1\rvert<7.3\times10^{-3}$ and $\sigma_t/\bar t\in[0.993,1.003]$; Sr$^+$ slow rate memory $s<0.27$ with the tested persistence model; $^{60}$Co dead time $t_g<1.3\times10^4\,\mathrm s$ or $5.3\times10^{-5}$ of mean lifetime. | Restricts nonexponential waiting-time shape, static or persistent heterogeneity, and a universal age-triggered minimum wait only in the stated models. The core exponential has $\sigma_t/\bar t=1$. |
| **DERIVED**, 45 | Quantum-dot trace clock shift $[-6,+10]\,\mathrm{ps}$; blocked-emission gate $t_g<53\,\mathrm{ps}$, with the assumed $5\,\mathrm{ps}$ systematic included. | Restricts those two specified onset models for that emitter. It does not give a threshold in energy units. The blocked-gate preference was only about $1.4\sigma$. |
| **DERIVED**, 46 | $^{57}$Fe forward-scattering clock shift $[-0.53,+0.42]\,\mathrm{ns}$, fractional ceiling $3.7\times10^{-3}$; blocked release $t_g<1.9\,\mathrm{ns}$, fractional ceiling $1.4\times10^{-2}$. | Bounds clock shift and hard gate relative to the synchrotron prompt in that nuclear response. Collective forward scattering is not a direct observation of an isolated particle unlocking. |
| **DERIVED**, 31/47 | FIRAS $\lvert\eta\rvert<9.0\times10^{-5}$ for a frequency-independent equilibrium bias at $x=\hbar\omega/k_BT\sim1$–$11$; Stage 31's $\sim50$ ppm rms is a spectral residual, not a pointwise bound on arbitrary $\eta(f)$. Johnson-noise thermometry $5.9\times10^{-14}$ for the Stage 47 low-$x$, constant-bias response projected through the particular readout; Raman/diode comparison $5.3\times10^{-6}$ under the repo's conservative interpretation. An $\eta=\varepsilon x$ variant reads as a temperature rescaling to first order, with $\lvert\varepsilon\rvert<6.2\times10^{-6}$ in the repo's cross-method comparison. | These bounds concern a bias in equilibrium forward/backward rates, not driven Schmitt hysteresis. The Johnson value also relies on the repo's $1/f$-like perturbation and weighting assumption; it is not a general frequency-local or universal bound on driven loops. |
| **FITTED**, 48/49 | CO ladder fit with temperature freed gives $\lvert\eta\rvert<2.3\times10^{-5}$ per rung under its convention; contextual closure $5.6\times10^{-5}$ is automatically satisfied by a common population model, arbitrary level dependence is constrained only in path sums at about $\pm2\times10^{-4}$. Transmon thermometry gives $\lvert\eta\rvert<6.6\times10^{-2}$ at $x\approx1.2$–$4.3$. | The apparent CO $4.3\sigma$ trend with temperature held fixed falls to $1.2\sigma$ when freed. The transmon result does not supersede FIRAS in its regime. |

**CITED** Grundmann et al. (2020), *Science* **370**, pp. 339–341, Eqs. (2)–(3), Fig. 3 and note 19, provides the phase-delay conversion, light-travel baseline, histogram methodology, and absolute-delay limitation. The public histogram is [Zenodo 3899920](https://doi.org/10.5281/zenodo.3899920). **DERIVED** $R/v_{\rm ph}=9204.3\,\mathrm{zs}$ per unit fringe coordinate at the Stage 52 inputs, where $v_{\rm ph}=E_e/p_e$. $R/c=246.8\,\mathrm{zs}$ for $R=0.74$ Å. The $504.5\,\mathrm{zs}$ fit used Gaussian plus background and a wider window. Its factor $2.04$ relative to $R/c$ is coincidentally close to both the phase/group velocity ratio and the span/slope factor. Likewise $268.7\,\mathrm{zs}$ from a different frequency convention is numerically near the $\sim270\,\mathrm{zs}$ nondipole comparison. Neither proximity establishes a common cause.

**CITED** Additional experimental locations used by the repo: Marx et al. (2013), *Physical Review Letters* **110**, 254801, abstract; Bollinger et al. (1989), *Physical Review Letters* **63**, p. 1031; Ossiander et al. (2017), *Nature Physics* **13**, pp. 280–285, Table 1; Norman et al. (1988), *Physical Review Letters* **60**, p. 2246, Table II; Ollivier et al. (2020), *ACS Photonics* **7**, p. 1050, Fig. 2b; Marx-Glowna et al. (2021), *Journal of Synchrotron Radiation* **28**, p. 120, Fig. 6; Qu et al. (2017), *Metrologia* **54**, p. 549; Li et al. (2025), *Science Advances* **11**, eadz6560, supplementary Table S1. **UNVERIFIED** Where exact equation or page and a transfer to this model are missing here, use the stage source CSVs and treat the corresponding historical inference as a repo result, not as a new published limit.


**DERIVED** For the separate equilibrium variant, define $x=\hbar\omega/(k_BT)$ and let $\eta$ multiply the thermal lock-to-unlock ratio after degeneracy factors are accounted for. The radiation occupation from the Stage 31 mapping is $n=[(1+\eta)e^x-1]^{-1}$, equivalent to a dimensionless chemical potential $\mu=\ln(1+\eta)$ in the exponent. For $x,|\eta|\ll1$, the leading expression is $n\simeq1/(x+\eta)$, which explains the low-$x$ leverage and also requires a physically admissible occupation over the domain of the model. A universal multiplicative bias around closed paths is only consistent with equilibrium if the logarithmic biases sum to zero around every cycle; applying the same bias naively to all upward transitions fails a triangle containing both two-rung and matching direct transitions. The Stage 48 ladder result assumes a rung-by-rung potential, not an unrestricted violation of detailed balance.

### Frequency constraints on the optional force-noise readout

**ASSUMED** $\lambda$ below denotes the rate in the repo's CSL-equivalent convention, with the detector transfer and mass geometry of each stage. A mechanical mode constrains a narrow weighted band near its mode, a broadband spectrum constrains resolved bins, and bulk heating and counting experiments constrain integrals. References and numerical inputs for individual rows are in the matching [stage source ledgers](../hyperfield/results/) and [Stage 54 aggregate](../hyperfield/results/stage54.json). The appendix enumerates all 52 entries in that aggregate, plus the later-assembled LPF envelope as a separate curve because Stage 54's aggregate omits it. These are reported nominal constraints, with heterogeneous confidence conventions and some indicative estimates. They are not a simultaneous 95% confidence contour.

**DERIVED** Main exclusions and exceptions:

- **DERIVED, Stages 31/33/53** LPF 2016 gives $2.96\times10^{-8}\,\mathrm{s^{-1}}$ at $0.7$–$20\,\mathrm{mHz}$; the Stage 33 conversion has its minimum $3.31\times10^{-9}\,\mathrm{s^{-1}}$ around $2$–$10\,\mathrm{mHz}$. Stage 53 uses all 17 clean L4 differential-acceleration segments and reaches $1.28\times10^{-8}\,\mathrm{s^{-1}}$ at $31.7$–$35.6\,\mathrm{mHz}$, then weakens to about $1.9\times10^{-3}\,\mathrm{s^{-1}}$ at the last accepted band, $1.002$–$1.125\,\mathrm{Hz}$. Stage 53's reported cut is approximately $1.12\,\mathrm{Hz}$. The former $0.19$–$9.3\,\mathrm{Hz}$ gap is covered only up to that cut; $1.12$–$9.3\,\mathrm{Hz}$ remains open in this map. Above about $50\,\mathrm{mHz}$ the bound rises roughly as $f^4$ from displacement readout noise; the $1\,\mathrm{Hz}$ line and $70\,\mathrm{mHz}$ bump are retained in conservative ceilings, not treated as candidates. The lowest per-segment $\lambda$ envelope uses marginal PSD bounds, with no simultaneous coverage claim.
- **DERIVED, Stages 33/34/54** The Rb expansion reaches the $\sim0.19\,\mathrm{Hz}$ kernel half-power point, not a sharp cutoff. The Zheng levitated-sphere validation gives $4.16\times10^{-7}\,\mathrm{s^{-1}}$ against the paper's $10^{-6.4}\,\mathrm{s^{-1}}$ (ratio $1.04$); the older $6.7\times10^{-7}$ review convention is conservative. LIGO's O3 minimum is $2.16\times10^{-6}\,\mathrm{s^{-1}}$ at $35.2\,\mathrm{Hz}$, deteriorating rapidly at high frequency; the analysis extends to $5\,\mathrm{kHz}$ but quantified calibration is stated for $10$–$2000\,\mathrm{Hz}$. The AURIGA point, ions and LIGO's upper band are weak.
- **DERIVED, Stage 37** At $8174\,\mathrm{Hz}$ Vinante's reported force floor was $1.87\,\mathrm{aN^2/Hz}$, a ceiling near $1.9$–$2.0\times10^{-8}\,\mathrm{s^{-1}}$ if entirely assigned to CSL-equivalent noise. At $3532.7\,\mathrm{Hz}$ the later multilayer instrument limits it to $2.0\times10^{-10}\,\mathrm{s^{-1}}$. For a line centered at the old mode and taking $\lambda_*=10^{-7.7}\,\mathrm{s^{-1}}$, substituting into the Lorentzian yields $b\le |8174-3532.7|\sqrt{r/(1-r)}\approx467\,\mathrm{Hz}$ with $r=2.0\times10^{-10}/10^{-7.7}$. The Gaussian counterpart has $\sigma\lesssim1530\,\mathrm{Hz}$. **ASSUMED** The suggested $b\gtrsim10\,\mathrm{Hz}$ is from taking the apparently flat $\sim20\,\mathrm{Hz}$ excess window literally, not a published line-width measurement. **FITTED** Turning the pulse tube on increases the floor by $0.71\pm0.26\,\mathrm{aN^2/Hz}$ ($2.8\sigma$ statistical); the later data are $1.9\sigma$ below the same unscaled noise floor. An environmental contribution is plausible and a white CSL origin at the excess strength fails other controls.
- **DERIVED, Stages 36/40–43/54** Stage 54 reports 25 narrow lines across its nominal gap convention, with only $0.019$ of $2.75$ frequency decades covered (about $0.7\%$ at $\pm0.1\%$ bookkeeping width). The first new points below $121\,\mathrm{kHz}$ are $25.22\,\mathrm{kHz}$ at $3.6\times10^{-6}$ and $100.5\,\mathrm{kHz}$ at $5.6\times10^{-6}\,\mathrm{s^{-1}}$, far above $10^{-7}$. The main empty stretches are $8.18$–$25.2\,\mathrm{kHz}$, $305$–$599\,\mathrm{kHz}$, and $1.80$–$4.60\,\mathrm{MHz}$ apart from the weak limiting ion line. Interior points at or below $10^{-7}$ are $242\,\mathrm{kHz}$, $1.406$, $1.486$ and $1.8\,\mathrm{MHz}$, with the $1.406\,\mathrm{MHz}$ point bath-implied. The $8174\,\mathrm{Hz}$ ceiling is on the lower edge. Seis and Gisler agree within about a factor $1.6$, compatible with similar geometries and thermal budgets. Above $4.6\,\mathrm{MHz}$ the $10.56$, $15.1$, $15.6$, and $65.4\,\mathrm{MHz}$ points and the $5.0486\,\mathrm{GHz}$ HBAR point do not constitute continuous coverage. Stage 43 raises the reliability of the transfer with bulk plus surface heating; the corrected quartz points remain bath-implied.
- **DERIVED, Stages 33–35** Bulk heating has a response centered near $9\,\mathrm{GHz}$ and excludes a broadband white rate above the repo's approximately $3.1\times10^{-11}\,\mathrm{s^{-1}}$ equivalent under its heating budget. The tail bound rises from about $4.0\times10^{-10}\,\mathrm{s^{-1}}$ near $13\,\mathrm{GHz}$ to roughly unity near $36.6\,\mathrm{GHz}$. Below roughly $12.4\,\mathrm{keV}$, neutral-atom coherent cancellation suppresses radiation for this mass-proportional long-smearing model; weak diffuse free-charge CMB/X-ray limits do not replace laboratory sensitivity. Neither the cryogenic low-energy excess nor the previously suggested neutral calorimeter resolves this window in this model. IGEX, Majorana and underground Ge constrain keV–MeV response; $100\,\mathrm{keV}$–$1\,\mathrm{MeV}$ remains without the proposed public differential count conversion. The non-mass-proportional Majorana result is a different model. Its electron-only and coherent-nucleus variants are not independent likelihoods from the same exposure.
- **UNVERIFIED, Stages 33/38/50/54** Neptune, intergalactic gas heating, neutron-star heating, CMB-proton heating, H nondissociation, and SQUID supercurrent limits were collated without adequate spectral response to impose on a narrow line. Earlier nuclear-spread scoping was superseded in part by Stage 44. N$_2$/CO intensity asymmetries are not timing data. In-band candidates dropped by Stage 54 lacked calibrated mode temperatures, measured heating or adequate geometry: Tao, Beccari, the $11.6\,\mathrm{kHz}$ Fogliano device, and Héritier's lower-frequency modes. Stage 41's archive-access obstruction was superseded by Stage 53's public-archive retrieval.

## 3. What survives jointly

**DERIVED (conditional on §1)** Let $B=425\,\mathrm{zs}$ use the matched fit, or $489\,\mathrm{zs}$ use the conservative envelope. With the gate active during the modeled H$_2$ events, both forks require $0\le Q_-<Q_+$, $\tau_0>|d|/2$, and, in the linear odd approximation,
$$
|d|<B[1+(\omega_{\rm e}\tau_0)^2].
$$
**DERIVED (fork A)** Add the visibility ceiling $\tau_0\lesssim0.76$–$0.90\,\mathrm{as}$, conditional on the observed visibility and baseline above. The corresponding linearized $d$ ceilings at the two visibility benchmarks are approximately $730$–$850\,\mathrm{zs}$ matched or $840$–$980\,\mathrm{zs}$ conservative, from $B/V_{\rm obs}$. At larger $|dx|$, use the exact phase and visibility factors from §1 and fit the same angular bins before treating those projections as exclusions.

**DERIVED (fork B)** The shared-phase model has no visibility ceiling on $\tau_0$. At $\tau_0=\hbar/E_e$, the leading matched ceiling doubles to $850\,\mathrm{zs}$, or $978\,\mathrm{zs}$ with the conservative limit. The $850\,\mathrm{zs}$ example is conditional on that illustrative $\tau_0$, not a global ceiling. Larger $\tau_0$ can further weaken the timing bound, subject to positivity and the exact phase law. Both forks admit $d=0$. The gate loads have no numerical calibration in the repo. The single-photon nondemolition result still excludes a universal high deposited-energy hurdle, while a driven gate with intact undriven detailed balance remains possible.

**DERIVED** For the spectral extension, let $W_i(f)\ge0$ denote the response converting spectral rate to measured excess in experiment $i$, and let $B_i>0$ be that experiment's excess budget in observable units. With no canceled noise contributions the surviving region, conditional on these response models, is
$$
f_0>0,\quad b>0,\quad
0\le\lambda_*\le
\inf_{i:\,I_i>0}\frac{B_i}{I_i(f_0,b)},\qquad
I_i(f_0,b)=\int_0^\infty W_i(f)\,
\frac{\bar g_i(f)\,df}{1+((f-f_0)/b)^2}.
$$
For a narrow calibrated mode centered at $f_i$, the repo approximation reduces this to $\lambda_*\bar g_i\le\lambda^{\max}_i[1+((f_i-f_0)/b)^2]$. For a resolved spectrum use this formula per bin or convolve with its bin response. For an integrated calorimeter or photon count use the integral, rather than taking the quoted white-noise limit as a pointwise exclusion of the same size. A gap means no local test there; a broad line centered in the gap can still be excluded by its wings.

**DERIVED (weak-coupling limit only, conditional on §1)** For either assumed H$_2$ phase-to-release mapping, an allowed core point can trivially be combined with sufficiently small positive spectral strength under the finite roster of positive budgets and finite calibrated responses at fixed $f_0,b$. This is continuity toward zero coupling and supplies no evidence for the hypothesis. The displayed upper envelope is positive whenever at least one experiment responds, and the full-exposure choice is conservative for $0\le\bar g_i\le1$. A sufficient bound is $\lambda_*<\min_i B_i/\int W_i\,df$ for the finite same-coupling set, since the Lorentzian and exposure are each at most unity. This argument requires no choice of a favorable resonance center or fit to an anomaly. It does not provide a lower bound, discovery probability, or a fundamental seven-parameter solution shared by atoms and solids. The repository quotes scalar equivalents and selected response curves; a joint statistically calibrated contour and gate calibration are **UNVERIFIED**. The open region shrinks or may disappear if a microscopic model fixes a minimum visible coupling, ties $\tau$ to its resonance, or asserts a specific unavoidable bias. No such relation is derived in the repo.

## 4. Three tests with a stated failure condition

| Status and experiment | Signal and target | Existing data and what a null would kill |
|---|---|---|
| **ASSUMED** Multi-energy H$_2$ coincidence photoionization, with event-level momentum and bond length; primary new setting $1600\,\mathrm{eV}$ paired with the $800\,\mathrm{eV}$ run | **DERIVED** A $100\,\mathrm{zs}$ planning benchmark gives an odd residual $100\cos\beta\,\mathrm{zs}$, a $200\,\mathrm{zs}$ end-to-end span, and about $0.0109\cos\beta$ in the Stage 52 fringe-coordinate units at the original energy. Preregister a common bond-length intercept at each energy, with the cross-energy dependence fixed by $d_1(E_e)=d/[1+(E_e\tau_0/\hbar)^2]$, versus a full nondipole/correlation baseline, and require $\sigma(d_1)\lesssim20\,\mathrm{zs}$ for nominal $5\sigma$ at the benchmark, including systematics. Megahertz or equilibrium noise parameters are irrelevant to this fit. | **DERIVED** The public histogram is a diagnostic, not decisive. Stage 32 forecasts $169\,\mathrm{zs}$ for a single-energy shape discriminator and $85\,\mathrm{zs}$ with the $1600\,\mathrm{eV}$ setting under its assumptions. Reducing $85\,\mathrm{zs}$ to $20\,\mathrm{zs}$ requires about $(85/20)^2\simeq18$ times the Stage 32 two-energy forecast's effective statistics if statistics dominate. **ASSUMED** Allow roughly $50$–$100$ times the original beamtime as a planning estimate accounting for the lower cross-section at $1600\,\mathrm{eV}$; the actual accepted coincidence yield and run duration are **UNVERIFIED**. **FITTED** Stage 52 had a $\pm32\,\mathrm{zs}$ fit-model systematic. More events alone cannot reach the $20\,\mathrm{zs}$ total target unless that uncertainty is reduced. **DERIVED** At the corrected $171\,\mathrm{zs}$ residual those imply expected sensitivities of about $1.0\sigma$ and $2.0\sigma$. Ask for the event file and fit settings, but obtain independent runs and acceptance reversals. A null at the design precision kills a stable $100\,\mathrm{zs}$ odd contribution in this channel under the preregistered mapping. |
| **ASSUMED** Retunable cryogenic mass-loaded cantilever or diamond nanoladder with calibrated mode thermometry, two mass profiles, and randomized cryocooler states | **DERIVED** Sweep across $8174\,\mathrm{Hz}$ before jumping to the $8$–$25\,\mathrm{kHz}$ gap; use steps smaller than the candidate $10\,\mathrm{Hz}$ HWHM in the initial sweep. The historical peak target is $1.87\,\mathrm{aN^2/Hz}$ or roughly $2.0\times10^{-8}\,\mathrm{s^{-1}}$ in the original geometry. **ASSUMED** A useful noise precision goal is $\sim0.19\,\mathrm{aN^2/Hz}$ in equivalent original geometry. **DERIVED** At the Stage 54 nanoladder conversion, reaching $10^{-7}\,\mathrm{s^{-1}}$ requires $\Gamma\lesssim6.4\times10^3\,\mathrm{s^{-1}}$, with $Q/T\gtrsim2\times10^7\,\mathrm{K^{-1}}$ for thermal-limited readout. | **DERIVED** Raw 2017/2020 spectra can test line flatness, vibration correlation and background, but the two fixed frequencies cannot establish an independently reproducible spectrum. A null with sensitivity below the size of the old excess across the candidate line excludes that continuous Lorentzian interpretation for the assumed mass coupling. It cannot eliminate arbitrarily weak or off-band lines. |
| **ASSUMED** Synchronized single-emitter gate test: quantum-dot excitation/readout with a calibrated Stark or strain load, two histories terminating at the same plateau, plus matched no-gate quantum evolution controls | **DERIVED** If $Q_-<q<Q_+$, approaching the plateau from above versus below predicts release over an observation window $T$ of $1-e^{-T/\tau_\theta}$ versus zero for the gated fraction. At $T=\tau_\theta$ the ideal difference is about $0.632$ times that fraction. Toggle pulse direction and ramp speed, and measure the prompt and emitted waveform on the same clock. **ASSUMED** Use an extra $10\,\mathrm{ps}$ clock shift as a distinct design benchmark and aim for $\lesssim2\,\mathrm{ps}$ total error. The gate test has no predicted signal size: the gated fraction, threshold loads, and conversion to a clock shift are unknown. The $10\,\mathrm{ps}$ value is only a design benchmark. | **DERIVED** Stage 45's digitized single trace constrains a shift or enforced dead time without comparing histories. The trace cannot measure both thresholds by itself. A matched-history null excludes a calibrated active threshold window and any extra shift at the stated sensitivity; it leaves a zeptosecond-scale or inactive channel open. Quantum-dot memory effects, including reservoir, detuning and pulse-history effects, are confounders and must be quantified before a loop can be attributed to the gate. |

## 5. Prior work and scope of novelty

- **CITED** Wigner (1955), *Physical Review* **98**, pp. 145–147, discusses the energy derivative of scattering phase as a time observable. Grundmann et al. (2020), *Science* **370**, pp. 339–341, Eqs. (2)–(3), measures angle-dependent relative phase timing from light travel across H$_2$. Ossiander et al. (2017), *Nature Physics* **13**, pp. 280–285, Table 1, resolves channel-dependent photoemission timing. **DERIVED** A smooth angle-dependent phase delay is therefore not itself novel. The proposed residual is a site-odd term with a specified bond-length and drive/history test after the ordinary phase model is controlled.
- **CITED** Minev et al. (2019), *Nature* **570**, pp. 200–204, [arXiv:1803.00545v3](https://arxiv.org/abs/1803.00545), Figs. 3–4, shows that a continuously monitored quantum jump can be tracked and reversed, consistent with conditional quantum trajectories. **DERIVED** This is the closest experimental precedent for the gradual-hook picture. Continuous conditioned evolution and genuine modification of unconditional measurement statistics are different tests. A slow trace or recovery under feedback alone does not establish a new mechanism.
- **CITED** Ghirardi, Rimini and Weber (1986), *Physical Review D* **34**, pp. 470–491, formulates stochastic localization; Ghirardi, Pearle and Rimini (1990), *Physical Review A* **42**, pp. 78–89, develops continuous spontaneous localization for identical particles. Carlesso, Ferialdi and Bassi (2018), *European Physical Journal D* **72**, article 159, sections on noninterferometric bounds, analyzes colored-collapse spectra. **DERIVED** Colored CSL is the closest precedent for the optional force-noise calculation. A fitted Lorentzian rate in mass-proportional CSL units is an implementation of that kind of probe. The hook-specific extra content would have to be an independently derived relationship to angle, threshold and release timing.
- **CITED** Diósi (1989), *Physical Review A* **40**, pp. 1165–1174, and Penrose (1996), *General Relativity and Gravitation* **28**, pp. 581–600, propose gravity-associated state reduction. Penrose (1996), p. 581 (abstract), gives the finite lifetime scale $\tau_{\rm DP}\sim\hbar/E_G$, with $E_G$ the gravitational self-energy of the difference between the mass distributions. **DERIVED** The present core supplies neither a gravitational self-energy functional nor a prediction of the Newton coupling. Converting the repository's separate DP/CQ limits to an angle-dependent phase or microscopic release time would require additional dynamics.
- **CITED** Zurek (2003), *Reviews of Modern Physics* **75**, pp. 715–775, reviews environment-induced decoherence and pointer states. **DERIVED** Tracing out a bath can suppress accessible coherences and yield history-dependent responses while leaving the enlarged dynamics ordinary. Driven bistability, delayed transitions, lifetime spread, and classical hysteresis are not individually distinguishing. **UNVERIFIED** Priority for the complete proposed angle × gate × release relation is not established because a microscopic relation has not been specified and a comprehensive priority search was outside these comparisons.

- **CITED** Vos et al. (2018), *Science* **360**, p. 1326, measures orientation-dependent stereo Wigner delay in CO. Nandi et al. (2020), *Science Advances* **6**, eaba7762, studies molecular shape-resonance emission timing. Holzmeier et al. (2021), *Nature Communications* **12**, 7343, studies how shape resonances affect the angular dependence of molecular photoionization delays.
- **CITED** Schultze et al. (2010), *Science* **328**, p. 1658, and Klünder et al. (2011), *Physical Review Letters* **106**, 143002, report attosecond photoemission timing. Pazourek, Nagele and Burgdörfer (2015), *Reviews of Modern Physics* **87**, p. 765, reviews the interpretation of those observables. Alongside Wigner's phase derivative above, Smith (1960), *Physical Review* **118**, p. 349, gives the scattering lifetime matrix.
- **CITED** Eckle et al. (2008), *Science* **322**, p. 1525, Landsman et al. (2014), *Optica* **1**, p. 343, and Sainadh et al. (2019), *Nature* **568**, p. 75, use angular streaking to study tunneling time. Ramos et al. (2020), *Nature* **583**, p. 529, uses a Larmor clock to address time within a barrier. **DERIVED** These timing observables use different clocks and do not validate the release-phase identification in §1.
- **CITED** Pearle (1976), *Physical Review D* **13**, p. 857, models dynamical state-vector reduction. Adler and Bassi (2007), *Journal of Physics A* **40**, pp. 15083–15098, develop collapse dynamics with non-white noise. Bassi et al. (2013), *Reviews of Modern Physics* **85**, p. 471, reviews dynamical collapse and experimental tests. Gammaitoni et al. (1998), *Reviews of Modern Physics* **70**, p. 223, reviews stochastic resonance, including noise-assisted response to periodic forcing.
- **DERIVED** The distinct test proposed here is a signed, tagged-site timing residual that follows the specified electron-energy dependence and either carries the independent-release visibility penalty or follows the shared-phase prediction. Linking that residual to a calibrated loading and unloading history would add a second discriminator. Finite duration, angular delays, resonant behavior and dynamical collapse all have precedents. **UNVERIFIED** A microscopic law tying the timing, visibility and gate into a single quantitative prediction, and priority for such a law, remain open.

## 6. Strongest objection

**DERIVED** The strongest objection is identifiability. An experiment measures an outgoing wave amplitude, a conditional trajectory, or excess device noise. A positive residual phase can be reparameterized as an ordinary scattering phase; a threshold loop as driven open-system memory; and a narrow force-noise peak as environmental noise. A strictly common stochastic delay cancels from relative phase; an independent release can reduce visibility even when the mean delays coincide. An additive theory with freely vanishing strength and unknown thresholds can agree with any finite list of upper limits. Describing a localized wave packet as a particle and its propagation as release does not by itself change quantum mechanics.

**DERIVED** The core is useful as a restricted falsifiable *phenomenology* only where its inputs and predictions are fixed before data are opened. A credible next step is to demonstrate the same signed residual with equivalent-site reversals, reject conventional $R$ and energy scaling in independent data, and test a history-conditioned gate with independently calibrated load. For a unified account, deduce one shared parameter relation that predicts another observable without fitting it, specify where the energy and momentum of any added stochastic kicks go, and supply a normalized dynamics compatible with known symmetries and equilibrium when undriven. **UNVERIFIED** No such derivation is presently in the sourced record, and the observations do not compel an extra channel.


## 7. Power-law release-noise channel

**ASSUMED** Use a one-sided spectral rate in the appendix's mass-proportional CSL-equivalent normalization, with the same reference smearing and $\bar g=1$:

$$
S(f)=A\left(\frac{f}{f_{\rm ref}}\right)^{-\alpha},\qquad f_{\rm ref}=1\,\mathrm{Hz}.
$$

Here $f$ is frequency in Hz, $A\ge0$ is the rate in $\mathrm{s^{-1}}$ at the reference frequency, and $\alpha$ is a dimensionless slope. $S$ here is neither the measured displacement power spectral density nor a spatial fractal dimension of the H₂ histogram. The numerical screen below is restricted to the listed experimental windows. **DERIVED** A nonzero power law over all positive frequencies has divergent total integrated power at one or both ends. A physical model needs breaks or cutoffs. No cutoff is estimated from these data, and promoting cutoffs to adjustable parameters adds further freedom.

**FITTED (free, unestimated)** The core gives neither $A$ nor $\alpha$, and there is no independent measurement of a release-noise component to fix either before comparison. The calculations below screen the full family without optimizing either parameter or assigning a best fit. The $\alpha=1$ column is the definition of the $1/f$ benchmark, not a slope inferred from the data. **DERIVED** The core has four parameters and the power law adds two, for six with the power law alone, or nine if both the three-parameter Lorentzian and the power law are present. The bounds below allocate the entire budget to the power law. An added Lorentzian consumes some of that budget.

**CITED** Carlesso, Ferialdi and Bassi (2018), *European Physical Journal D* **72**, 159, Eq. (17), uses frequency-weighted bulk heating for colored CSL. **DERIVED** For an actual nonnegative response $W_i$ and measured upper budget $B_i$, the conditional allowed region is $0\le A\le \min_i B_i/K_i(\alpha)$, where

$$
K_i(\alpha)=\int W_i(f)(f/f_{\rm ref})^{-\alpha}\,df.
$$

For a band quoted as an equivalent constant rate $L_i=B_i/\int W_i\,df$, the upper edge is $L_i/\langle (f/f_{\rm ref})^{-\alpha}\rangle_{W_i}$. For a resolved mode at $f_i$, this reduces to $A\le L_i(f_i/f_{\rm ref})^\alpha$ under the narrow-mode approximation. Integrated X-ray counts, expansion and heating cannot be recast as pointwise limits using only the appendix. For a positive kernel known to be supported in $[l,h]$, dividing $L_i$ by the maximum spectral shape on that interval gives a sufficient ceiling and dividing by the minimum gives a necessary ceiling. The Rb half-power point and the bulk-heating half-maximum interval are sensitivity descriptors. They are not hard bounds on kernel support.

### Numerical comparison with the complete appendix

**ASSUMED** To put numbers on the plane with the published summaries, approximate each unresolved band by a normalized response uniform in $\ln f$, use mode centers for narrow mechanical lines, and use the resolved curve samples where available. This is a reproducible screen, and the resulting boundaries apply to this surrogate response. The Stage 53 line-ceiling envelope can also be screened with the same bin rule. The band limit in each bin is tighter and controls the combined envelope there, without treating the two estimators as independent. The Stage 35 tail is screened bin by bin using its own quoted budgets, not summed as independent measurements. The Stage 34 numerical edge uses the quantified calibration range. The bins above that range remain indicative. The appendix's non-mass-proportional Majorana row belongs to a different coupling and is excluded from the joint envelope. The electron-only Majorana row is screened as a weaker alternate projection of the same exposure rather than as an independent likelihood.

**DERIVED** For each log-uniform band the divisor is

$$
H_\alpha(l,h)=
\begin{cases}
1 & \alpha=0,\
$$2pt]
\displaystyle\frac{(l/f_{\rm ref})^{-\alpha}-(h/f_{\rm ref})^{-\alpha}}{\alpha\ln(h/l)} & \alpha\ne0.
\end{cases}
$$

Take $A_{\max}(\alpha)=\min_i L_i/H_\alpha$, using the point rule for modes or individual samples. This specifies the numerical allowed region for any slope in the stated surrogate. It is an intersection of quoted ceilings, with no joint confidence level assigned. The switches below were checked on $-3\le\alpha\le4$.

| Status | Slope interval, approximate | Active upper edge on $A$ |
|---|---:|---|
| **DERIVED** | $-3$ to $-1.24$ | Underground Ge, Donadi MeV band |
| **DERIVED** | $-1.24$ to $0.269$ | Majorana X-rays, mass-proportional coherent nuclei |
| **DERIVED** | $0.269$ to $0.626$ | Stage 33 interpolated LPF 2018 spectrum, with the active frequency changing across the sampled curve |
| **DERIVED** | $0.626$ to $4$ | Rb free expansion with the assumed log-uniform band and the repository's $10^{-6}\,\mathrm{Hz}$ plotting floor |

**DERIVED** Representative $A_{\max}$ values in $\mathrm{s^{-1}}$ at $1\,\mathrm{Hz}$ are $1.77\times10^{-54}$ at $\alpha=-2$, $4.19\times10^{-34}$ at $-1$, $4.94\times10^{-15}$ at $0$, $1.50\times10^{-10}$ at $0.5$, $6.20\times10^{-13}$ at $1$, and $1.24\times10^{-18}$ at $2$. **UNVERIFIED** The steep positive-slope edge is particularly sensitive to the assumed infrared response and plotting floor. It is not a measured Rb bound on $1/f$ noise. The exact experimental contour and exact locations of the edge switches require the true response kernels, exposure window, spectral cutoffs and covariance. The table does not establish those inputs.

**DERIVED** The following audit identifies every Stage 54 aggregate entry by its order in the machine-readable source and shows its individual $A$ ceiling at $\alpha=1$. The original appendix remains the source for frequencies and limits. Curve entries use the underlying samples or bins instead of applying their best limit uniformly across the band.

| Status | Entry | Probe | $A$ ceiling at $\alpha=1$ ($\mathrm{s^{-1}}$) |
|---|---:|---|---:|
| **DERIVED** | 1 | LISA Pathfinder 2016 (Helou 2017) | $7.20\times10^{-11}$ |
| **DERIVED** | 2 | LISA Pathfinder 2018 (Armano 2018 ASD -> lambda, DERIVED) | $6.65\times10^{-12}$ (curve or bins) |
| **DERIVED** | 3 | Cold atoms, 87Rb pK free expansion (Kovachy data, Bilardello 2016) | $6.20\times10^{-13}$ |
| **DERIVED** | 4 | Diamagnetic levitated micro-sphere 9.3 Hz (Zheng 2020) | $6.23\times10^{-6}$ |
| **DERIVED** | 5 | Diamagnetic levitated micro-sphere 12.9 Hz (Zheng 2020) | $8.64\times10^{-6}$ |
| **DERIVED** | 6 | LIGO test-mass force noise 30-35 Hz | $3.24\times10^{-4}$ |
| **DERIVED** | 7 | AURIGA bar 931 Hz | $4.28\times10^{1}$ |
| **DERIVED** | 8 | Nanocantilever 3084 Hz (Vinante 2016) | $1.54\times10^{-4}$ |
| **DERIVED** | 9 | Multilayer cantilever 3532.7 Hz (Vinante 2020) | $7.07\times10^{-7}$ |
| **DERIVED** | 10 | Cantilever 8174 Hz (Vinante 2017) | $1.55\times10^{-4}$ |
| **DERIVED** | 11 | 40Ca+ heating 4.6 MHz (Kaufmann 2017) | $6.28\times10^{6}$ |
| **DERIVED** | 12 | 88Sr+ heating 0.60 MHz (Bruzewicz 2015) | $4.23\times10^{5}$ |
| **DERIVED** | 13 | 88Sr+ heating 0.80 MHz (Bruzewicz 2015) | $6.16\times10^{5}$ |
| **DERIVED** | 14 | 88Sr+ heating 1.00 MHz (Bruzewicz 2015) | $5.79\times10^{5}$ |
| **DERIVED** | 15 | 88Sr+ heating 1.30 MHz (Bruzewicz 2015) | $6.58\times10^{5}$ |
| **DERIVED** | 16 | 88Sr+ heating 1.50 MHz (Bruzewicz 2015) | $9.73\times10^{5}$ |
| **DERIVED** | 17 | Bulk heating of solids (10 pW/kg residual; Adler-Vinante 2018) | $2.55\times10^{-1}$ |
| **DERIVED** | 18 | HBAR phonon mode 5.961 GHz (Schrinski 2023) | $1.42\times10^{6}$ |
| **DERIVED** | 19 | IGEX Ge X-rays 14.5-48.5 keV (Piscicchia 2017) | $4.11\times10^{7}$ |
| **DERIVED** | 20 | Majorana 19-100 keV, m-p electrons only | $1.64\times10^{6}$ |
| **DERIVED** | 21 | Majorana 19-100 keV, m-p coherent nuclei | $4.65\times10^{4}$ |
| **DERIVED** | 22 | Majorana 19-100 keV, NON-mass-proportional | Excluded from mass-proportional envelope |
| **DERIVED** | 23 | Ge underground 1000-3800 keV (Donadi 2021) | $2.28\times10^{8}$ |
| **DERIVED** | 24 | LIGO O3a L1 / O3b H1 / O1 GW150914 force noise (DERIVED) | $1.00\times10^{-4}$ (curve or bins) |
| **DERIVED** | 25 | Bulk-heating tail above 13 GHz, band-limited (0.05-dec bins, lambda<=1 shown) | $5.55\times10^{0}$ (curve or bins) |
| **DERIVED** | 26 | CXB total 1-2 keV vs free-charge emission of diffuse ionized gas (Hickox 2006) | $1.83\times10^{19}$ |
| **DERIVED** | 27 | CXB total 2-8 keV vs free-charge emission of diffuse ionized gas (Hickox 2006) | $3.29\times10^{19}$ |
| **DERIVED** | 28 | COBE/FIRAS 2-21 cm^-1 residual vs free-charge emission (Fixsen 1996) | $4.65\times10^{19}$ |
| **DERIVED** | 29 | HBAR mode 404, 5.0486 GHz (Omahen 2025) [bulk+surface] | $7.82\times10^{3}$ |
| **DERIVED** | 30 | SiN membrane 1.14 MHz, 20 nm (Rossi 2018), heating n_th Gamma_m at T~10 K | $1.36\times10^{0}$ |
| **DERIVED** | 31 | SiN membrane 1.486 MHz, 63 nm + Al pad (Seis 2022), 1/t_coh = n_th Gamma_m | $4.30\times10^{-3}$ |
| **DERIVED** | 32 | SiN membrane 242 kHz, Device I (1,1) (Yuan 2015) | $4.27\times10^{-3}$ |
| **DERIVED** | 33 | SiN membrane 121 kHz, Device II (1,1) (Yuan 2015) | $2.33\times10^{-2}$ |
| **DERIVED** | 34 | SiN membrane 193 kHz, Device II (1,2) (Yuan 2015) | $5.92\times10^{-2}$ |
| **DERIVED** | 35 | SiN membrane 192 kHz, Device II (2,1) (Yuan 2015) | $2.44\times10^{0}$ |
| **DERIVED** | 36 | Al drum 10.56 MHz (Teufel 2011) | $1.46\times10^{2}$ |
| **DERIVED** | 37 | Al vacuum-gap drum 1.8 MHz (Youssefi 2023) | $6.57\times10^{-2}$ |
| **DERIVED** | 38 | Al drum 15.1 MHz, 500 uK (Cattiaux 2021) | $2.98\times10^{1}$ |
| **DERIVED** | 39 | quartz BAW OT5 15.6 MHz, 20 mK (Goryachev 2012) [bulk+surface] | $7.34\times10^{1}$ |
| **DERIVED** | 40 | quartz BAW OT21 65.4 MHz, 20 mK (Goryachev 2012) [bulk+surface] | $1.79\times10^{3}$ |
| **DERIVED** | 41 | HBAR ~6.65 GHz (Chu 2017) [bulk+surface] | $8.88\times10^{6}$ |
| **DERIVED** | 42 | Diamond nanoladder cantilever 25.22 kHz (Heritier 2018) | $9.12\times10^{-2}$ |
| **DERIVED** | 43 | SiN trampoline 41.4 kHz, 44 nm (Reinhardt 2016) | $9.99\times10^{-1}$ |
| **DERIVED** | 44 | Levitated silica nanosphere z 77.6 kHz, cryogenic (Tebbenjohanns 2021) | $1.72\times10^{0}$ |
| **DERIVED** | 45 | Levitated silica nanosphere z 104 kHz, room-T (Magrini 2021) | $3.97\times10^{0}$ |
| **DERIVED** | 46 | Levitated silica nanosphere 100.5 kHz, R = 68.0 nm (Jain 2016) | $5.63\times10^{-1}$ |
| **DERIVED** | 47 | Levitated silica nanosphere 148.8 kHz, R = 52.7 nm (Jain 2016) | $2.63\times10^{0}$ |
| **DERIVED** | 48 | Levitated silica nanosphere 151.2 kHz, R = 71.6 nm (Jain 2016) | $1.71\times10^{0}$ |
| **DERIVED** | 49 | Levitated silica nanosphere 158.8 kHz, R = 68.0 nm (Jain 2016) | $2.77\times10^{0}$ |
| **DERIVED** | 50 | Levitated silica nanosphere x 305 kHz (Delic 2020) | $4.39\times10^{1}$ |
| **DERIVED** | 51 | SiN tethered membrane 140 kHz, 20 nm, room T (Norte 2016) | $7.65\times10^{0}$ |
| **DERIVED** | 52 | Soft-clamped SiN string 1.406 MHz, 20 nm, 46 mK (Gisler 2022) | $6.54\times10^{-3}$ |

**DERIVED** In addition to those 52 entries, both estimators in every one of the 31 accepted Stage 53 LPF bins were evaluated with $L_i/H_1$. The strongest band ceiling is about $4.29\times10^{-10}\,\mathrm{s^{-1}}$ at $0.0317$–$0.0356\,\mathrm{Hz}$; the corresponding line-ceiling screen is $4.65\times10^{-10}\,\mathrm{s^{-1}}$. The six above-cut bins were screened separately at the same benchmark. Their best band ceiling is about $3.89\times10^{-3}\,\mathrm{s^{-1}}$, and they stay outside the robust envelope. Without the Rb surrogate the $\alpha=1$ screen is set by the interpolated Stage 33 LPF curve at $6.65\times10^{-12}\,\mathrm{s^{-1}}$. The stronger Stage 33 low-frequency limits retain their interpolated status.

**DERIVED** A preregistered power law with fixed nonzero $A$ and fixed $\alpha$ is easier to test across widely separated bands because the same slope and normalization apply everywhere. Giving the hypothesis independent freedom in $A$ and $\alpha$, especially an arbitrarily small $A$ or adjustable breaks, makes the hypothesis as a whole harder to falsify. No agreement between a chosen slope and a spectral feature counts as independent evidence. **UNVERIFIED** Neither a release-noise slope nor its physical link to a spatial fractal is established.

## 8. Nearby mass at the H₂ detector

**DERIVED — predicted delay shift per kilogram at 10 cm: undetermined.** The release model contains no coupling from an external source mass $M$, gravitational potential, or tidal field into $\tau_0$, $d$ or $Q_\pm$. In the leading timing law, defining $F=[1+(\omega_{\rm e}\tau_0)^2]^{-1}$ and holding electron energy fixed gives

$$
\frac{\partial\Delta t}{\partial M}=\cos\beta\left[
F\,\frac{\partial d}{\partial M}
-2d\,\omega_{\rm e}^{2}\tau_0 F^2\,\frac{\partial\tau_0}{\partial M}\right].
$$

Neither derivative is supplied. The written model has zero explicit response to $M$ at fixed parameters. **UNVERIFIED** A nonzero nearby-mass extension has no predicted size in $\mathrm{zs/kg}$ and is not yet falsifiable. A mass-swapping null could instead bound a separately defined empirical slope.

**ASSUMED** Use the same empty outer housing, the same housing with an aluminum insert, and a housing with a dense pure tungsten insert. Specify and magnetically screen the actual tungsten stock and mount. Avoid ferromagnetic tungsten alloys or binders. Match outer shape, position and grounding. Use the measured center-to-interaction distance, with $r=0.10\,\mathrm m$ as the calculation benchmark, and replace the point approximation with the measured mass geometry when needed.

**DERIVED** The Newtonian comparison for a compact body is $\Phi/M=-G/r\simeq-6.67\times10^{-10}\,\mathrm{m^2\,s^{-2}\,kg^{-1}}$. The potential itself is common to the two ends of the molecule to leading order. Across $R=0.74$ Å aligned with the field, the potential difference per kilogram is approximately $GR/r^2=4.94\times10^{-19}\,\mathrm{m^2\,s^{-2}\,kg^{-1}}$. A weak-field clock-rate comparison gives a maximum fractional differential scale of about $5.5\times10^{-36}\,\mathrm{kg^{-1}}$. Multiplying by the pre-existing $896$ zs reference time gives an illustrative scale of $4.9\times10^{-33}\,\mathrm{zs/kg}$. This is a clock-rate scale comparison, not a calculation of the spectrometer or of a Newtonian shift in the measured birth delay. A constant potential alone cannot supply a locally measured differential delay. **UNVERIFIED** A full ordinary-physics fringe prediction including the apparatus geometry and gravitational transport has not been calculated here.

**FITTED** Stage 52 gives $\pm35.8$ zs statistical and $\pm32$ zs fit-model uncertainty. **DERIVED** Their quadrature is about $48.0$ zs if treated as independent, some 34 orders above the illustrative differential clock scale for a $1$ kg swap. This comparison supplies no predicted signal scale for a new mass coupling and no detection claim. **ASSUMED** If a mass coupling is proposed, specify its sign, distance law and minimum effect before opening the swap data. Do not select a distance or a slope because it coincides with a fringe feature.

| Status | Confounder | Control |
|---|---|---|
| **ASSUMED** | Magnetic fields or gradients affecting electron trajectories | Use screened tungsten, aluminum and nonmagnetic mounts. Map the field for each insert with the spectrometer energized and track a fixed calibration feature. |
| **ASSUMED** | Electrostatic fields, charging and patch potentials | Use the identical grounded conductive exterior, check ground continuity, and repeat the trajectory calibration after each swap. |
| **ASSUMED** | Thermal drift from different heat capacity or conduction | Temperature-match inserts before a run, log the chamber and mount temperatures, and use the same stabilization criterion for all three cases. |
| **ASSUMED** | Vibration, mechanical loading and alignment | Support the insert independently of the instrument, monitor acceleration and alignment, and include blinded empty-housing sham swaps. |
| **ASSUMED** | Blackbody radiation from changed temperature, emissivity or view factor | Use the same temperature-monitored exterior and radiation shield for every insert, preserve the line of sight and view factors, and control temperature before acquisition. |

**ASSUMED** Randomize and blind swaps, repeat them on both sides, and run calibration checks without changing detector settings. A repeatable mass-correlated residual would motivate testing its side and distance dependence after bounding the confounders. **UNVERIFIED** No nearby-mass measurement or release-specific calibration has been made in this report.

## Appendix: numerical frequency inventory

**DERIVED** The following inventory transcribes the Stage 54 aggregate at readable precision, leaving the originals as the authoritative machine-readable values. Where a curve is named, the number is the best value in the listed band, rather than a uniform limit. The mechanical lines use the repo's $\pm0.1\%$ plotting widths. The stage labels point to the experiment-specific derivations and paper locations. The non-mass-proportional alternative and indicative estimates carry their own assumptions.

| Status | Stage | Device or coupling | Hz interval or mode center | Repo equivalent $\lambda_{\max}$ ($\mathrm{s^{-1}}$) | Note |
|---|---:|---|---:|---:|---|
| **DERIVED** (repo inventory) | 31 | LISA Pathfinder 2016 (Helou 2017) | $7\times10^{-4}$–$2\times10^{-2}$ | $2.96\times10^{-8}$ |  |
| **DERIVED** (repo inventory) | 33 | LISA Pathfinder 2018 (Armano 2018 ASD -> lambda, DERIVED) | $2\times10^{-5}$–$3\times10^{-2}$ | $3.31\times10^{-9}$ | min 3.31e-09 at 2-10 mHz; 3.94e-06 at 20 uHz; review quotes 3.8e-9 |
| **DERIVED** (repo inventory) | 33 | Cold atoms, 87Rb pK free expansion (Kovachy data, Bilardello 2016) | $1\times10^{-6}$–$1.89\times10^{-1}$ | $5.1\times10^{-8}$ | responds to all f below the half-power point; lower edge = plot edge |
| **DERIVED** (repo inventory) | 33 | Diamagnetic levitated micro-sphere 9.3 Hz (Zheng 2020) | $9.3$ | $6.7\times10^{-7}$ |  |
| **DERIVED** (repo inventory) | 33 | Diamagnetic levitated micro-sphere 12.9 Hz (Zheng 2020) | $12.9$ | $6.7\times10^{-7}$ |  |
| **DERIVED** (repo inventory) | 33 | LIGO test-mass force noise 30-35 Hz | $3\times10^{1}$–$3.5\times10^{1}$ | $1\times10^{-5}$ |  |
| **DERIVED** (repo inventory) | 33 | AURIGA bar 931 Hz | $931$ | $4.6\times10^{-2}$ |  |
| **DERIVED** (repo inventory) | 33 | Nanocantilever 3084 Hz (Vinante 2016) | $3084$ | $5\times10^{-8}$ |  |
| **DERIVED** (repo inventory) | 33 | Multilayer cantilever 3532.7 Hz (Vinante 2020) | $3532.7$ | $2\times10^{-10}$ |  |
| **DERIVED** (repo inventory) | 31 | Cantilever 8174 Hz (Vinante 2017) | $8174.01$ | $1.9\times10^{-8}$ |  |
| **DERIVED** (repo inventory) | 33 | 40Ca+ heating 4.6 MHz (Kaufmann 2017) | $4600000$ | $1.37$ |  |
| **DERIVED** (repo inventory) | 33 | 88Sr+ heating 0.60 MHz (Bruzewicz 2015) | $600000$ | $7.04\times10^{-1}$ |  |
| **DERIVED** (repo inventory) | 33 | 88Sr+ heating 0.80 MHz (Bruzewicz 2015) | $800000$ | $7.7\times10^{-1}$ |  |
| **DERIVED** (repo inventory) | 33 | 88Sr+ heating 1.00 MHz (Bruzewicz 2015) | $1000000$ | $5.79\times10^{-1}$ |  |
| **DERIVED** (repo inventory) | 33 | 88Sr+ heating 1.30 MHz (Bruzewicz 2015) | $1300000$ | $5.06\times10^{-1}$ |  |
| **DERIVED** (repo inventory) | 33 | 88Sr+ heating 1.50 MHz (Bruzewicz 2015) | $1500000$ | $6.48\times10^{-1}$ |  |
| **DERIVED** (repo inventory) | 33 | Bulk heating of solids (10 pW/kg residual; Adler-Vinante 2018) | $5.55\times10^{9}$–$1.3\times10^{10}$ | $3.1\times10^{-11}$ | peak 9.00e+09 Hz; review quotes 3.3e-11 |
| **DERIVED** (repo inventory) | 33 | HBAR phonon mode 5.961 GHz (Schrinski 2023) | $5961000000$ | $2.39\times10^{-4}$ |  |
| **DERIVED** (repo inventory) | 33 | IGEX Ge X-rays 14.5-48.5 keV (Piscicchia 2017) | $3.51\times10^{18}$–$1.17\times10^{19}$ | $6.8\times10^{-12}$ |  |
| **DERIVED** (repo inventory) | 31 | Majorana 19-100 keV, m-p electrons only | $4.59\times10^{18}$–$2.42\times10^{19}$ | $1.74\times10^{-13}$ |  |
| **DERIVED** (repo inventory) | 31 | Majorana 19-100 keV, m-p coherent nuclei | $4.59\times10^{18}$–$2.42\times10^{19}$ | $4.94\times10^{-15}$ |  |
| **DERIVED** (repo inventory) | 31 | Majorana 19-100 keV, NON-mass-proportional | $4.59\times10^{18}$–$2.42\times10^{19}$ | $5.15\times10^{-20}$ | different model (n-m-p); not comparable with mechanical bounds; excluded from coverage stats |
| **DERIVED** (repo inventory) | 33 | Ge underground 1000-3800 keV (Donadi 2021) | $2.42\times10^{20}$–$9.19\times10^{20}$ | $5.2\times10^{-13}$ |  |
| **DERIVED** (repo inventory) | 34 | LIGO O3a L1 / O3b H1 / O1 GW150914 force noise (DERIVED) | $1\times10^{1}$–$5\times10^{3}$ | $2.16\times10^{-6}$ | [Response curve](../hyperfield/results/stage34_ligo_curve.csv) |
| **DERIVED** (repo inventory) | 35 | Bulk-heating tail above 13 GHz, band-limited (0.05-dec bins, lambda<=1 shown) | $1.3\times10^{10}$–$3.66\times10^{10}$ | $4.04\times10^{-10}$ | [Response curve](../hyperfield/results/stage35_bulk_tail_curve.csv) |
| **DERIVED** (repo inventory) | 35 | CXB total 1-2 keV vs free-charge emission of diffuse ionized gas (Hickox 2006) | $2.42\times10^{17}$–$4.84\times10^{17}$ | $5.47\times10^{1}$ | all measured CXB attributed (+2 sigma); free charges spaced >> lambda_k -> incoherent |
| **DERIVED** (repo inventory) | 35 | CXB total 2-8 keV vs free-charge emission of diffuse ionized gas (Hickox 2006) | $4.84\times10^{17}$–$1.93\times10^{18}$ | $3.68\times10^{1}$ | all measured CXB attributed (+2 sigma); free charges spaced >> lambda_k -> incoherent |
| **DERIVED** (repo inventory) | 35 | COBE/FIRAS 2-21 cm^-1 residual vs free-charge emission (Fixsen 1996) | $6\times10^{10}$–$6.3\times10^{11}$ | $2.98\times10^{8}$ | weak: CMB is bright |
| **DERIVED** (repo inventory) | 43 | HBAR mode 404, 5.0486 GHz (Omahen 2025) [bulk+surface] | $5048630000$ | $1.55\times10^{-6}$ | paper own convention: 5.7e-8; with paper 1-sigma P_max: 7.16e-07 |
| **DERIVED** (repo inventory) | 36 | SiN membrane 1.14 MHz, 20 nm (Rossi 2018), heating n_th Gamma_m at T~10 K | $1140000$ | $1.19\times10^{-6}$ | T only "~10 K" in text: factor-2 uncertainty |
| **DERIVED** (repo inventory) | 36 | SiN membrane 1.486 MHz, 63 nm + Al pad (Seis 2022), 1/t_coh = n_th Gamma_m | $1486000$ | $2.89\times10^{-9}$ |  |
| **DERIVED** (repo inventory) | 40 | SiN membrane 242 kHz, Device I (1,1) (Yuan 2015) | $242000$ | $1.77\times10^{-8}$ | Q CITED (14 mK) |
| **DERIVED** (repo inventory) | 40 | SiN membrane 121 kHz, Device II (1,1) (Yuan 2015) | $121000$ | $1.92\times10^{-7}$ | Q >= 5e6 CITED (lower bound used -> conservative) |
| **DERIVED** (repo inventory) | 40 | SiN membrane 193 kHz, Device II (1,2) (Yuan 2015) | $193000$ | $3.07\times10^{-7}$ | Q >= 5e6 CITED (lower bound used -> conservative) |
| **DERIVED** (repo inventory) | 40 | SiN membrane 192 kHz, Device II (2,1) (Yuan 2015) | $192000$ | $1.27\times10^{-5}$ | Q = 1.2e5 CITED; mode T of Device II ASSUMED to apply |
| **DERIVED** (repo inventory) | 41 | Al drum 10.56 MHz (Teufel 2011) | $10560000$ | $1.39\times10^{-5}$ | first membrane point in the 4.6 MHz-5.6 GHz window |
| **DERIVED** (repo inventory) | 41 | Al vacuum-gap drum 1.8 MHz (Youssefi 2023) | $1800000$ | $3.65\times10^{-8}$ | thickness from FEM caption (ASSUMED equal to device) |
| **DERIVED** (repo inventory) | 42 | Al drum 15.1 MHz, 500 uK (Cattiaux 2021) | $15100000$ | $1.97\times10^{-6}$ | measured occupancy; rho_Al from Teufel |
| **DERIVED** (repo inventory) | 43 | quartz BAW OT5 15.6 MHz, 20 mK (Goryachev 2012) [bulk+surface] | $15600000$ | $4.71\times10^{-6}$ | bath-implied heating n_TH*w/Q; no occupancy measurement |
| **DERIVED** (repo inventory) | 43 | quartz BAW OT21 65.4 MHz, 20 mK (Goryachev 2012) [bulk+surface] | $65400000$ | $2.74\times10^{-5}$ | bath-implied heating n_TH*w/Q; no occupancy measurement |
| **DERIVED** (repo inventory) | 43 | HBAR ~6.65 GHz (Chu 2017) [bulk+surface] | $6650000000$ | $1.34\times10^{-3}$ | above Omahen 5.05 GHz; weak |
| **DERIVED** (repo inventory) | 54 | Diamond nanoladder cantilever 25.22 kHz (Heritier 2018) | $25220$ | $3.62\times10^{-6}$ | tier A; heat from T,Q_lo = 1.4e+05; from F_th(+1sig), m(-1sig) = 2.3e+05 s^-1; rails treated as independent beams (conservative) |
| **DERIVED** (repo inventory) | 54 | SiN trampoline 41.4 kHz, 44 nm (Reinhardt 2016) | $41400$ | $2.41\times10^{-5}$ | tier A; Q = pi f tau = 3.71e+07; measured x_rms 165+-5 pm vs thermal 161 pm and measured S_F spectrum within 5% of thermal |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere z 77.6 kHz, cryogenic (Tebbenjohanns 2021) | $77600$ | $2.22\times10^{-5}$ | tier A; m used 1e-18 kg, R 50.5 nm |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere z 104 kHz, room-T (Magrini 2021) | $104000$ | $3.82\times10^{-5}$ | tier A; no uncertainties quoted on Gamma_ba, Gamma_th |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere 100.5 kHz, R = 68.0 nm (Jain 2016) | $100500$ | $5.6\times10^{-6}$ | tier A; Gamma quoted in s^-1 (no 2pi), cf. their Eq. (11) 13.0 kHz |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere 148.8 kHz, R = 52.7 nm (Jain 2016) | $148800$ | $1.77\times10^{-5}$ | tier A; Gamma quoted in s^-1 (no 2pi), cf. their Eq. (11) 13.0 kHz |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere 151.2 kHz, R = 71.6 nm (Jain 2016) | $151200$ | $1.13\times10^{-5}$ | tier A; Gamma quoted in s^-1 (no 2pi), cf. their Eq. (11) 13.0 kHz |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere 158.8 kHz, R = 68.0 nm (Jain 2016) | $158800$ | $1.74\times10^{-5}$ | tier A; Gamma quoted in s^-1 (no 2pi), cf. their Eq. (11) 13.0 kHz |
| **DERIVED** (repo inventory) | 54 | Levitated silica nanosphere x 305 kHz (Delic 2020) | $305000$ | $1.44\times10^{-4}$ | tier A;  |
| **DERIVED** (repo inventory) | 54 | SiN tethered membrane 140 kHz, 20 nm, room T (Norte 2016) | $140000$ | $5.46\times10^{-5}$ | tier B; tether width not in fetched text: plate limit ASSUMED |
| **DERIVED** (repo inventory) | 54 | Soft-clamped SiN string 1.406 MHz, 20 nm, 46 mK (Gisler 2022) | $1406000$ | $4.65\times10^{-9}$ | tier B; check: kB T/(hbar Q) = 2.62 s^-1 |

**DERIVED** Stage 53 LPF L4 envelope from the [complete CSV](../hyperfield/results/stage53_lpf_lambda_vs_f.csv). The marginal 95% band and line columns use the repo's separate estimators. The [above-cut values](../hyperfield/results/stage53_lpf_lambda_indicative_above_cut.csv) are indicative and excluded from the robust envelope.

| Status | Band (Hz) | Band rate ($\mathrm{s^{-1}}$) | Line rate ($\mathrm{s^{-1}}$) |
|---|---:|---:|---:|
| **DERIVED** | $0.031698$–$0.035566$ | $1.28\times10^{-8}$ | $1.39\times10^{-8}$ |
| **DERIVED** | $0.035566$–$0.039905$ | $1.59\times10^{-8}$ | $1.83\times10^{-8}$ |
| **DERIVED** | $0.039905$–$0.044774$ | $2\times10^{-8}$ | $2.48\times10^{-8}$ |
| **DERIVED** | $0.044774$–$0.050238$ | $2.85\times10^{-8}$ | $3.6\times10^{-8}$ |
| **DERIVED** | $0.050238$–$0.056368$ | $3.91\times10^{-8}$ | $4.68\times10^{-8}$ |
| **DERIVED** | $0.056368$–$0.063246$ | $5.55\times10^{-8}$ | $7.28\times10^{-8}$ |
| **DERIVED** | $0.063246$–$0.070963$ | $1.3\times10^{-7}$ | $2.79\times10^{-7}$ |
| **DERIVED** | $0.070963$–$0.079621$ | $1.31\times10^{-7}$ | $1.68\times10^{-7}$ |
| **DERIVED** | $0.079621$–$0.089337$ | $1.76\times10^{-7}$ | $2.25\times10^{-7}$ |
| **DERIVED** | $0.089337$–$0.10024$ | $2.51\times10^{-7}$ | $3.02\times10^{-7}$ |
| **DERIVED** | $0.10024$–$0.11247$ | $3.9\times10^{-7}$ | $5.07\times10^{-7}$ |
| **DERIVED** | $0.11247$–$0.12619$ | $5.73\times10^{-7}$ | $7.58\times10^{-7}$ |
| **DERIVED** | $0.12619$–$0.14159$ | $9.19\times10^{-7}$ | $1.23\times10^{-6}$ |
| **DERIVED** | $0.14159$–$0.15887$ | $1.44\times10^{-6}$ | $1.88\times10^{-6}$ |
| **DERIVED** | $0.15887$–$0.17825$ | $2.23\times10^{-6}$ | $2.92\times10^{-6}$ |
| **DERIVED** | $0.17825$–$0.2$ | $3.53\times10^{-6}$ | $4.58\times10^{-6}$ |
| **DERIVED** | $0.2$–$0.2244$ | $5.63\times10^{-6}$ | $7.42\times10^{-6}$ |
| **DERIVED** | $0.2244$–$0.25179$ | $8.91\times10^{-6}$ | $1.18\times10^{-5}$ |
| **DERIVED** | $0.25179$–$0.28251$ | $1.42\times10^{-5}$ | $1.99\times10^{-5}$ |
| **DERIVED** | $0.28251$–$0.31698$ | $2.24\times10^{-5}$ | $2.85\times10^{-5}$ |
| **DERIVED** | $0.31698$–$0.35566$ | $3.45\times10^{-5}$ | $4.54\times10^{-5}$ |
| **DERIVED** | $0.35566$–$0.39905$ | $5.41\times10^{-5}$ | $7.27\times10^{-5}$ |
| **DERIVED** | $0.39905$–$0.44774$ | $8.25\times10^{-5}$ | $1.11\times10^{-4}$ |
| **DERIVED** | $0.44774$–$0.50238$ | $1.27\times10^{-4}$ | $1.7\times10^{-4}$ |
| **DERIVED** | $0.50238$–$0.56368$ | $1.93\times10^{-4}$ | $2.52\times10^{-4}$ |
| **DERIVED** | $0.56368$–$0.63246$ | $2.92\times10^{-4}$ | $3.9\times10^{-4}$ |
| **DERIVED** | $0.63246$–$0.70963$ | $4.32\times10^{-4}$ | $5.62\times10^{-4}$ |
| **DERIVED** | $0.70963$–$0.79621$ | $6.4\times10^{-4}$ | $8.25\times10^{-4}$ |
| **DERIVED** | $0.79621$–$0.89337$ | $9.28\times10^{-4}$ | $1.28\times10^{-3}$ |
| **DERIVED** | $0.89337$–$1.0024$ | $1.37\times10^{-3}$ | $3.24\times10^{-3}$ |
| **DERIVED** | $1.0024$–$1.1247$ | $1.9\times10^{-3}$ | $2.43\times10^{-3}$ |

**DERIVED** Reproduction checks for the forecast numbers above: $100/9204.3=0.01086$, $171/169=1.01$, $171/85=2.01$, $1-e^{-1}=0.63212$, and $(6.4\times10^{10})(10^{-7})=6.4\times10^3$. The cantilever two-frequency ratio gives a Lorentzian half-width of about $467\,\mathrm{Hz}$; the approximate $10\,\mathrm{Hz}$ lower scale hinges on the source's presumed flat excess window. Neither was fitted here.

**DERIVED** A final coincidence flag from Stage 51: the forward/backward visibility ratio $0.872$ lies near the counts ratio $0.873$. The ratios have no established causal relationship. Likewise no frequency is picked by overlap with Adler's proposed strength scale.
