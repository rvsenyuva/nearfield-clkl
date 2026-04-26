# Covariance-Domain Near-Field Channel Estimation under Hybrid Compression

**CL-KL: Curvature-Learning via KL Covariance Fitting**

MATLAB simulation code for the paper:

> R. V. Şenyuva, "Covariance-Domain Near-Field Channel Estimation under Hybrid
> Compression: USW/Fresnel Model, Curvature Learning, and KL Covariance Fitting,"
> *IEEE Transactions on Communications*, submitted April 2026.

---

## Overview

This repository contains the complete simulation codebase for the CL-KL near-field
channel estimator and all six comparison baselines. The code reproduces all 11
figures and tables in the paper from scratch using a fixed random seed.

**Key result:** CL-KL achieves the lowest channel NMSE among all six evaluated
methods at SNR ∈ {−5, 0, +5, +10} dB, operating exclusively on the
N_RF × N_RF compressed covariance (64 complex values at N_RF = 8) — a 64×
reduction in data volume versus full-array methods.

---

## Requirements

- MATLAB R2023b or later (tested: R2025b)
- Parallel Computing Toolbox (optional; required only for `parfor` mode)
- No additional toolboxes required

Tested on: Windows 11, AMD Ryzen 5 7500F (6-core), MATLAB R2025b (R2023b compatible).

---

## Repository Structure

```
nearfield_sim/
├── setup.m                    % Add all subfolders to MATLAB path (run first)
├── run_all.m                  % Master script: runs all 9 figure scripts
├── run_ablation.m             % Ablation study (4 configs × 3 SNR × N_MC=50)
├── test_modules.m             % Smoke tests for all estimators (10 tests)
├── test_fixes.m               % Regression tests for bug fixes (13 tests)
│
├── core/                      % Shared utilities
│   ├── nf_params.m            % Default parameter struct P
│   ├── nf_update_derived_pub.m% Recompute derived fields after changing P.M
│   ├── nf_gen_channel.m       % Generate random near-field channel
│   ├── nf_hybrid_combiner.m   % Random phase-only combining matrix W
│   ├── pregen.m               % Pre-generate N_MC trials (enables parfor)
│   ├── nf_fresnel_steer.m     % Fresnel steering vector a(theta, u)
│   ├── nf_fresnel_steer_du.m  % Derivative da/du for gradient computation
│   ├── nf_usw_steer.m         % Exact USW steering vector (truth model)
│   ├── save_results.m         % Save workspace results struct
│   ├── save_results_csv.m     % Append results to nf_simulation_results.csv
│   ├── aggregate.m            % Aggregate per-trial arrays to mean/std
│   ├── plot_fig.m             % Standard NMSE errorbar figure
│   ├── nf_add_legend.m        % Unified below-axes legend for multi-panel figs
│   ├── nf_export_fig.m        % Export figures to PDF (IEEE column widths)
│   └── setup.m                % (symlink — use root setup.m)
│
├── estimators/                % Channel estimation algorithms
│   ├── nf_clkl.m              % CL-KL (proposed): KL covariance fitting
│   ├── nf_psomp.m             % P-SOMP: beam-depth range sampling
│   ├── nf_zhang.m             % DL-OMP: two-subarray dictionary learning
│   ├── nf_music_tri.m         % MUSIC+Tri: subarray MUSIC + LS triangulation
│   ├── nf_dfrft_nomp.m        % DFrFT-NOMP: fractional Fourier transform NOMP
│   └── nf_bfsomp.m            % BF-SOMP: beam-focused SOMP (Hussain TWC 2025)
│
├── figures/                   % One script per paper figure
│   ├── run_fig2_nmse_snr.m    % Fig. 2: NMSE vs SNR (primary result)
│   ├── run_fig3_nmse_nrf.m    % Fig. 3: NMSE vs N_RF
│   ├── run_fig4_nmse_n.m      % Fig. 4: NMSE vs N (snapshots)
│   ├── run_fig5_nearfar.m     % Fig. 5: Near-to-far transition sweep
│   ├── run_fig6_runtime.m     % Fig. 6: Runtime vs M
│   ├── run_fig7_robustness.m  % Fig. 7: Fresnel model-mismatch robustness
│   ├── run_fig7b_source_robustness.m  % Fig. 7b: QPSK vs Gaussian pilots
│   ├── run_fig8_vard.m        % Fig. 8: NMSE vs number of paths d
│   └── run_fig9_convergence.m % Fig. 9: CL-KL convergence diagnostic
│
└── metrics/                   % Performance evaluation
    ├── nf_crb.m               % Compressed-domain Cramér–Rao Bound
    └── nf_metrics.m           % NMSE, RMSE(θ), RMSE(r), failure rate
```

---

## Quickstart

```matlab
% 1. Set path (run once per MATLAB session)
cd('/path/to/nearfield_sim')
setup

% 2. Verify installation
test_modules   % should print: 11 PASSED, 0 FAILED
test_fixes     % should print: 13 PASSED, 0 FAILED

% 3. Run all figures (parallel, ~5 minutes on 6-core machine)
run_all('parfor')

% 4. Run ablation study separately (~1 min with parfor)
run_ablation(true)
```

Output files written to the working directory:
- `nf_simulation_results.csv` — all figure results (30 columns)
- `nf_ablation_results.csv`   — ablation results (9 columns)
- `fig2_nmse_snr.pdf`, `fig3_nmse_nrf.pdf`, … — publication-ready figures

---

## Default Parameters

All figures use these defaults unless overridden (see `nf_params.m`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| `f_c` | 28 GHz | Carrier frequency |
| `M` | 64 | ULA elements (default) |
| `N_RF` | 8 | RF chains (default) |
| `N` | 64 | Snapshots per trial |
| `d` | 3 | Number of paths (default) |
| `theta` range | [20°, 60°] | Angle support |
| `r` range | [0.05, 1.0] × r_RD | Range support |
| `N_MC` | 400 | Monte Carlo trials (publication) |
| `SNR` sweep | −15 to +25 dB | 9-point sweep |
| Random seed | `rng(42,'twister')` | Fixed for reproducibility |

Rayleigh distance: r_RD = 2D²/λ = 21.26 m at M = 64, f_c = 28 GHz.

---

## Simulation Modes

```matlab
run_all             % Serial, full N_MC=400  (~20–40 min)
run_all('parfor')   % Parallel, full N_MC=400 (~5 min, 6 cores)
run_all('fast')     % Serial, N_MC=20 preview (~8 min)
run_all('fast parfor') % Parallel, N_MC=20 preview (~1 min)
```

Note: Fig. 6 (runtime measurement) always runs serially regardless of mode,
since parfor would distort wall-clock timing.

---

## Methods

| Method | Input | Data volume | Reference |
|--------|-------|-------------|-----------|
| **CL-KL** (proposed) | R̂_y | N_RF² = 64 values | This paper |
| P-SOMP | R̂_y | N_RF² = 64 values | Cui & Dai, TCOM 2022 + Hussain TWC 2025 |
| DL-OMP | X | M×N = 4096 values | Zhang et al., TCOM 2024 |
| MUSIC+Tri | X | M×N = 4096 values | Haghshenas et al., ICC 2025 |
| DFrFT-NOMP | X | M×N = 4096 values | Yang et al., WCL 2024 |
| BF-SOMP | X | M×N = 4096 values | Hussain et al., TWC 2025 |

**Filled markers** in all figures: compressed-domain input (CL-KL, P-SOMP).
**Open markers**: full-array input (64× more data than compressed methods).

---

## Reproducibility

All results in the paper are fully reproducible from this code:

```matlab
rng(42, 'twister')   % fixed seed — already set inside each run_fig*.m
run_all('parfor')    % regenerates all figures and nf_simulation_results.csv
run_ablation(true)   % regenerates ablation table (nf_ablation_results.csv)
```

Runtime on reference hardware (AMD Ryzen 5 7500F, 6 cores, 6000 MHz DDR5):
- Full suite (`run_all('parfor')`): ~5 minutes
- Ablation (`run_ablation(true)`): ~1 minute
- Test suites: ~30 seconds combined

---

## CSV Output Schema

### `nf_simulation_results.csv` (30 columns)

| Column | Description |
|--------|-------------|
| `timestamp` | ISO datetime of run |
| `figure_id` | Fig2, Fig3, …, Fig9 |
| `figure_name` | Human-readable figure name |
| `experiment_type` | Type of sweep |
| `fixed_M/N_RF/N/d/SNR_dB` | Fixed parameters |
| `sweep_variable/value` | Swept parameter |
| `truth_model` | exact_USW / gaussian / qpsk |
| `N_MC` | Trials used |
| `r_RD_m` | Rayleigh distance [m] |
| `method` | Algorithm name |
| `input_data` | R_hat or X |
| `NMSE_dB_mean/std` | Channel NMSE [dB] |
| `RMSE_theta_deg_mean` | Angle RMSE [deg] |
| `RMSE_r_m_mean` | Range RMSE [m] |
| `fail_rate_pct` | % trials with parameter error > threshold |
| `runtime_s` | Per-trial wall time [s] |
| `clkl_avg_iters` | Mean CL-KL iterations (CL-KL rows only) |
| `clkl_conv_pct` | % CL-KL converged (CL-KL rows only) |
| `clkl_N0_mean` | Mean N̂₀ estimate (CL-KL rows only) |
| `CRB_theta_deg` | Compressed-domain CRB angle [deg] |
| `CRB_r_m` | Compressed-domain CRB range [m] |
| `CRB_theta_full_deg` | Full-array CRB angle [deg] |
| `CRB_r_full_m` | Full-array CRB range [m] |
| `notes` | Free text |

### `nf_ablation_results.csv` (9 columns)

| Column | Description |
|--------|-------------|
| `timestamp` | ISO datetime |
| `config_name` | CL-KL (full) / w/o frozen N0 / w/o multi-start / w/o post-loop scan |
| `ablation_flags` | Which flag is active |
| `SNR_dB` | SNR level |
| `NMSE_dB_mean` | Channel NMSE [dB] |
| `RMSE_theta_deg_mean` | Angle RMSE [deg] |
| `RMSE_r_m_mean` | Range RMSE [m] |
| `fail_rate_pct` | Failure rate [%] |
| `N_MC` | Trials (always 50 for ablation) |

---

## Known Limitations

- **N_RF identifiability limit:** CL-KL and P-SOMP require d ≤ ⌊(N_RF−1)/2⌋
  (= 3 for N_RF = 8). NMSE degrades beyond this boundary.
- **High-SNR OMP warm-start:** At SNR = +25 dB, CL-KL's OMP residual deflation
  with a frozen-curvature approximation can select the wrong active set at
  extreme SNR, causing convergence to a local minimum. All baselines are
  unaffected. A Newton/Fisher-scoring refinement step is identified as
  future work.
- **Range accuracy:** CL-KL optimises channel NMSE, not individual range
  accuracy. The post-loop matched-filter scan partially recovers range
  accuracy but a gap to the compressed-domain CRB remains.

---

## Citation

```bibtex
@article{senyuva2026clkl,
  author  = {Şenyuva, Rıfat Volkan},
  title   = {Covariance-Domain Near-Field Channel Estimation under Hybrid
             Compression: {USW/Fresnel} Model, Curvature Learning, and
             {KL} Covariance Fitting},
  journal = {IEEE Transactions on Communications},
  year    = {2026},
  note    = {Submitted April 2026}
}
```

---

## License

Code released under the MIT License. See `LICENSE` for details.
