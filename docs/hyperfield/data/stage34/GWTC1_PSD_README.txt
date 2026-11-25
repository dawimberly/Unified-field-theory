The Power Spectral Densities (PSDs) reported in this data release have been used for all event specific analysis reported in GWTC-1.
Each event-specific file contains a column defining the frequency series (both bandwidth and resolution) used in the analysis of that event, as well as columns containing the PSD for each detector included in the analysis.
A description for how these PSDs were computed are given in Appendix B of the GWTC-1 paper, and the relevant section is quoted below:

	We estimate the power spectral density that enters the inner product using BayesWave [53, 54]. The PSD is modeled as a cubic spline for the broad-band structure and a sum of Lorentzians for the line features. A median PSD is computed from the resulting posterior PDF of PSDs, defined separately at each frequency. This PSD is expected to lead to more stable and reliable inference.

Each event has an associated set of PSDs that is valid only for the specific data segment that was included in the analysis, the PSDs are therefore not valid for analysis of any other data.
The data segments are defined as ending 2s after the times reported in Table I of the GWTC-1 paper, and have a total duration T defined as

	T = 1/(f[i] - f[i-1])

where f[i] represents the i-th value of the frequency array reported in the first column of each PSD data file.
The duration of the analysed data segment therefore sets the frequency resolution of the analysis, and subsequently also the frequency resolution of the PSD.
Using this frequecy resolution can make some features in the PSDs appear to be artificially truncated, but it is important to remember that this is indeed only an apparent truncation as all features are resolved appropriately for the analysis that's been undertaken.

Note that the LIGO and Virgo data are only properly calibrated in the frequency range between 10 Hz and 5 kHz [62, 63]; data at lower and higher frequencies should not be used for any astrophysical analysis. 
For the parameter estimation analysis reported in the GWTC-1 paper, only data in the frequency range from 20 Hz (30 Hz for GW170608 at Hanford, 23Hz for GW170817 and 16Hz for GW170818) and 1 kHz (2 kHz for GW170817) was used. 
Thus the PSDs in this data release only have values in that range.

These PSDs differ in apperance and behaviour from what is shown in the detector summary pages (https://www.gw-openscience.org/summary_pages/detector_status/).
The main difference is that the summary pages show Amplitude Spectral Densities (ASD = sqrt(PSD)) which are estimated using Welch's method (https://en.wikipedia.org/wiki/Welch%27s_method) from an average over a long duration of detector data.
The PSDs presented here are, as mentioned in the GWTC-1 paper quote above, estimated only from the data segment where the GW signal is believed to be observable using a model (consisting of Lorentzians for the line features and a cubic spline for the broadband features) implemented in BayesWave [53, 54].
This model includes the number of Lorenzians, as well as the behaviour of the spline, as free parameters in the analysis which generates a posterior probability distribution over the possible PSDs that are valid for whitening this specific data segment.
To capture the overall behaviour of this posterior distribution, the meadian PSD is computed for each frequency included in the analysis.

The PSDs for different events will differ significantly, this is primarily caused by the different behaviour of the detectors at different times, but also by the post-run subtraction of certain line features (calibration lines, harmonics of the power-grid frequency etc.) performed for the O2 data from the LIGO detectors.
This is described in Section II.C of the GWTC-1 paper, with the relevant section is quoted below:

	We subtracted several independent contributions to the instrumental noise from the data at both LIGO detectors [51]. For all of O2, the average increase in the BNS range from this noise subtraction process at LHO was ≈18% [51]. At LLO the noise subtraction process targeted narrow line features, resulting in a negligible increase in BNS range.
	Calibrated strain data from each interferometer was produced online for use in low-latency searches. Following the run, a final frequency-dependent calibration was generated for each interferometer. 
	For the LIGO instruments this final calibration benefitted from the use of post-run measurements and removal of instrumental lines. 

Again, to reiterate, the PSDs reported in this data release are only to be used for the specific data segments from which they were generated.

Overall, the model, configuration and use of the PSDs reported in this data release were chosen as they together fulfill the requirements of the Parameter Estimation analysis presented in the GWTC-1 paper (ie the likelihood that is defined assumes the whitened residuals, once a proposed gravitational wave model signal has been subtracted from the data, follows a normal distribution).


{The references below follow the same numbering as in the GWTC-1 paper}

[51] D. Davis, T. J. Massinger, A. P. Lundgren, J. C. Driggers, A. L. Urban, and L. K. Nuttall, “Improving the Sensitivity of Advanced LIGO Using Noise Subtraction,” (2018), arXiv:1809.05348 [astro-ph.IM]
[53] Neil J. Cornish and Tyson B. Littenberg, “BayesWave:Bayesian Inference for Gravitational Wave Bursts and Instrument Glitches,” Class. Quant. Grav. 32, 135012 (2015), arXiv:1410.3835 [gr-qc]
[54] Tyson B. Littenberg and Neil J. Cornish, “Bayesian inference for spectral estimation of gravitational wave detector noise,” Phys. Rev. D 91, 084034 (2015), arXiv:1410.3852 [gr-qc]
[62] Craig Cahillane et al. (LIGO Scientific Collaboration), “Calibration uncertainty for Advanced LIGO’s first and second observing runs,” Phys. Rev. D 96, 102001 (2017), arXiv:1708.03023 [astro-ph.IM]
[63] Aaron Viets et al., “Reconstructing the calibrated strain signal in the Advanced LIGO detectors,” Class. Quant. Grav. 35, 095015 (2018), arXiv:1710.09973 [astro-ph.IM] .
