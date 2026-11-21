# FFZ three-loop Newton matching review

This benchmark uses equal Yukawa couplings in a Wess–Zumino matter theory on an external metric. The fixed-\(\xi=1/6\) subtraction and the original dimensional continuation are part of the test.

## Start here

1. [Original conventions and two-loop calculation](ffz_two_loop_flat_space_checkpoint.md)
2. [Fixed-1/6 scheme conversion and the Newton moments](ffz_fixed_one_sixth_checkpoint.md)
3. [Graph inventory](ffz_newton_graphs_and_pole_mass_checkpoint.md) and [signed Weyl trace check](ffz_yukawa_spin_traces_check.py)
4. [Y improvement simple pole](ffz_three_loop_Y_improvement_checkpoint.md)
5. [Current scalar total and corrected QCc overlap](ffz_qcc_full_regulated_check.md), [QFf](ffz_qff_curvature_pole.md), and [latest CcFf regulated check](ffz_ccff_regulated_poles.md)

## Review target

The additive X curvature beta cancels through order h³ in the stated matter-source subtraction. The projected Newton test is open. The improvement-independent three-loop simple pole must equal 1/24 in units h²z/K to cancel the Y improvement at a_X=0, a_Y=−1. First verify cancellation of projected higher poles. The two curved FFff graphs and remaining counterterms are outstanding. The CcFf double pole is 1/144 and its simple-pole constant at L=0 is approximately −0.0247. The latest QCc checkpoint supersedes both earlier double-pole values and the earlier scalar sum. Use its corrected QCc double pole 1/72.

For Grok review, check the current QCc and CcFf subtractions, including the evanescent overlap and long-pair terms. Then check the signed FFff contractions before carrying them into the curvature calculation. Report exact results, numerical estimates, and unresolved terms separately. Preserve the distinction between a matter-source conversion and pure-metric matching.

## Reproduce

The checkpoint Python blocks run in their stated order. Requirements across the batch are Python 3, NumPy, SciPy, SymPy, and mpmath. From this directory run `python ffz_yukawa_spin_traces_check.py` for the explicit Wick and flat normalization check. For the regulated CcFf check, extract the Python block in the latest checkpoint and run it with `OPENBLAS_NUM_THREADS=1`. Archived checkpoints record the derivation history and can contain superseded numbers.

## Checkpoint archive

- [ffz_ccff_regulated_poles.md](ffz_ccff_regulated_poles.md)
- [ffz_ccff_long_pair.md](ffz_ccff_long_pair.md)
- [ffz_ccff_short_pair.md](ffz_ccff_short_pair.md)
- [ffz_ccff_projective_kernel.md](ffz_ccff_projective_kernel.md)
- [ffz_qff_curvature_pole.md](ffz_qff_curvature_pole.md)
- [ffz_qcc_full_regulated_check.md](ffz_qcc_full_regulated_check.md)
- [ffz_qcc_finite_scalar_sum.md](ffz_qcc_finite_scalar_sum.md)
- [ffz_qcc_missing_corner_checkpoint.md](ffz_qcc_missing_corner_checkpoint.md)
- [ffz_scalar_endpoint_residues.md](ffz_scalar_endpoint_residues.md)
- [ffz_scalar_projective_checkpoint.md](ffz_scalar_projective_checkpoint.md)
- [ffz_cubic_newton_checkpoint.md](ffz_cubic_newton_checkpoint.md)
- [ffz_quartic_newton_checkpoint.md](ffz_quartic_newton_checkpoint.md)
- [ffz_yukawa_spin_traces_check.py](ffz_yukawa_spin_traces_check.py)
- [ffz_newton_graphs_and_pole_mass_checkpoint.md](ffz_newton_graphs_and_pole_mass_checkpoint.md)
- [ffz_threshold_and_spinor_checkpoint.md](ffz_threshold_and_spinor_checkpoint.md)
- [ffz_subtracted_momentum_kernel.md](ffz_subtracted_momentum_kernel.md)
- [ffz_three_loop_Y_improvement_checkpoint.md](ffz_three_loop_Y_improvement_checkpoint.md)
- [ffz_newton_projection_q4_checkpoint.md](ffz_newton_projection_q4_checkpoint.md)
- [ffz_two_loop_momentum_slope.md](ffz_two_loop_momentum_slope.md)
- [ffz_fixed_one_sixth_checkpoint.md](ffz_fixed_one_sixth_checkpoint.md)
- [ffz_local_weyl_and_curvature_checkpoint.md](ffz_local_weyl_and_curvature_checkpoint.md)
- [ffz_two_loop_flat_space_checkpoint.md](ffz_two_loop_flat_space_checkpoint.md)
