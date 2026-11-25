# Composite-metric denominator control

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records composite-metric denominator control. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Calculation completed directly in the interactive session; not an automation request.

**DERIVED** Question: when does replacing an expectation of a composite ratio by a ratio of expectations lose terms? The CEV composite metric contains $P/V$. This control tests the approximation, not the actual renormalized CEV vacuum.

**DERIVED** Take N independent unit normal fields, $R=\sum (\phi _a^{2})/N, V=N(v+bR)$, with $v,b>0$. Let the numerator be independent of R with mean Np. Then the normalized ratio expectation is

$$
E[P/V]/p = \int _0^\infty \exp\left(-v t\right)(1+2bt/N)^{-N/2} dt.
$$

**DERIVED** The ratio-of-means approximation is $1/(v+b)$. Expanding the exact expectation at large N gives

$$
E[P/V]/p = 1/(v+b) + 2b^{2}/[N(v+b)^{3}] + O(N^{-2}).
$$

**DERIVED** Jensen's inequality fixes the correction to be positive in this positive ensemble. Twelve cases cover $N=16,64,256,1024$ and ($v,b)=(1$,.1),(1,1),(.2,1). Independent Simpson refinements at 4000,8000,16000 panels give maximum final refinement difference $4.494\times10^{-9}$. For $v=b=1$ the relative error of the leading correction decreases from 3.19% at $N=16$ to 0.0489% at $N=1024$.

**DERIVED** The integral tail is bounded by $\exp\left(-40\right)/v$ with the chosen cutoff. The calculation uses Python's standard library and an existing Simpson helper, included in the archive.

**DERIVED** At $v=0, b=1$, the exact inverse moment is $N/(N-2)$ for $N>2$ and diverges for $N\le 2$. Thus denominator-zero limits require an integrability check; an otherwise plausible large-N replacement is not uniformly justified across arbitrary singular branches.

**DERIVED** For a general correlated numerator and denominator, expanding about their means gives

$$
\begin{gathered}
E[P/V] = Pbar/Vbar - Cov(P,V)/Vbar^{2} \\
+ Pbar Var(V)/Vbar^{3} + higher joint cumulant terms.
\end{gathered}
$$

**DERIVED** These formulas identify connected insertions needed in a source-matched next-order composite-metric calculation. They do not determine their renormalized values or whether a physical scalar pole is shifted. Independence is a control assumption; an interacting field theory need not satisfy it. Positive Gaussian sampling does not implement the Lorentzian or dimensionally continued CEV measure, normal-ordering prescription, or constrained gauge quotient. No defect in the published leading-order calculation is demonstrated by this control.

**CITED** Published context: Batz, Erlich, Mrini, Background-Independent Composite Gravity, https://arxiv.org/pdf/2011.06541 (published DOI https://doi.org/10.1088/1361-6382/abefbd). Its composite metric and vacuum self-consistency motivate the question. The finite-dimensional integral and checks here are our control calculation, not a result attributed to that paper. — eq./page UNVERIFIED.

**ASSUMED** Next useful calculation: derive the corresponding connected composite insertions with the same regulator and operator prescription as the chosen field-theory vacuum. Do not substitute these positive-ensemble numbers into the CEV pole formula.
