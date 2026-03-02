# fdi_wglt

This folder implements a practical heterogeneous redundant IMU fault-detection pipeline
for a `2 INS + 1 ISIS` architecture (`n=9`, `m=3`) using whitened parity-space GLRT.

## Files

- `fdi_precompute_model.m`: offline precomputation (white transform, parity matrix,
  tolerable fault level, adaptive threshold).
- `fdi_generate_redundant_measurements.m`: builds synthetic redundant measurements
  from reference 3-axis signals and injects single or multi-window fault events.
- `fdi_run_sliding_wglt.m`: online sliding-window detection/isolation and fault-size estimate.

## Demo entry points

- `demos/test_FDI_WGLT_IRS_ISIS.m`: simple single-fault example.
- `demos/test_FDI_WGLT_fault_campaign.m`: batch campaign matching paper fault cases,
  with per-case visualization and metric summary.
