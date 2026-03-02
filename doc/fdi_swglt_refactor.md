# SWGLT FDI Refactor Notes

This repository now uses a corrected SWGLT pipeline for heterogeneous INS/ISIS
FDI processing.

## Core updates

- Replaced PPV/PCA denoising with sliding-window parity mean `P_bar`.
- Detection now uses energy statistic `FD_j = j * ||P_bar||^2`.
- Isolation now uses matched-filter statistics
  `FI_j(i) = j * (c_tilde_i' * P_bar)^2 / ||c_tilde_i||^2`.
- Added robust isolation gating:
  - margin ratio threshold `rho_threshold`
  - Bayesian posterior confidence `P_isol`
  - dwell confirmation count `N_dwell`
- Added configurable vote-based isolation confirmation `min_isolation_votes`
  to reduce overly conservative DNI-only outcomes when one confidence metric
  is temporarily weak.
- Adaptive thresholds are computed from non-central chi-square formulas using
  whitened matrices (`W`, `H_w`, `V_w`).

## New outputs from `fdi_swglt`

The algorithm returns a `results` struct with:

- `FD_energy`, `FD_max`
- `FI`
- `isolated`, `iso_status`, `f_hat`
- `posterior`
- `T_adaptive`, `T_adaptive_proj`
- `diagnostics`

## Integration changes

- `main_simulation.m` now consumes the SWGLT `results` struct.
- `evaluate_performance.m` supports DNI/CIR/MIR metrics for robust isolation.
- `plot_results.m` adds SWGLT margin and posterior diagnostic figures.

## Tuning demo

- `demos/test_SINS_FDI_SWGLT_tuning.m` provides an example grid sweep for
  `(rho_threshold, P_isol, N_dwell)` with `min_isolation_votes=2`.
