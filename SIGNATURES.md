# CL-KL Codebase — Function Signatures

Generated automatically. One entry per `.m` file.

> **Usage:** Run `setup.m` once per MATLAB session to add all subfolders to path.

## Core Utilities (`core/`)

### `aggregate.m`
```matlab
function [nmse_pt, nstd_pt, rth_pt, rr_pt, fail_pt, ...
```
AGGREGATE  Compute mean/std across Monte Carlo trials for all metrics. [nmse_pt, nstd_pt, rth_pt, rr_pt, fail_pt, ...

### `nf_add_legend.m`
```matlab
function lh = nf_add_legend(fig, methods, styles, colors, varargin)
```
NF_ADD_LEGEND  Place a single unified horizontal legend below all subplots. lh = nf_add_legend(fig, methods, styles, colors)

### `nf_export_fig.m`
```matlab
function nf_export_fig(fig_handle, filename, layout, varargin)
```
NF_EXPORT_FIG  Export a MATLAB figure to publication-quality PDF. nf_export_fig(fig_handle, filename, layout)

### `nf_fresnel_steer.m`
```matlab
function a = nf_fresnel_steer(theta, u, P)
```
NF_FRESNEL_STEER  Fresnel (second-order / chirp) near-field steering vector. a = nf_fresnel_steer(theta, u, P)

### `nf_fresnel_steer_du.m`
```matlab
function [da_du, c_i] = nf_fresnel_steer_du(theta, a, P)
```
NF_FRESNEL_STEER_DU  Derivative of Fresnel steering vector w.r.t. u = 1/r. [da_du, c_i] = nf_fresnel_steer_du(theta, a, P)

### `nf_gen_channel.m`
```matlab
function [X, H_true, S_true, theta_true, r_true, p_true, N0] = ...
```
NF_GEN_CHANNEL  Generate near-field multi-path channel snapshots. STEP 1 PHYSICAL REGIME AUDIT (TWC submission plan):

### `nf_hybrid_combiner.m`
```matlab
function [W, Y, R_hat] = nf_hybrid_combiner(X, P)
```
NF_HYBRID_COMBINER  Build random-phase combiner, compress, whiten, form sample covariance. [W, Y, R_hat] = nf_hybrid_combiner(X, P)

### `nf_params.m`
```matlab
function P = nf_params()
```
NF_PARAMS  Locked simulation defaults (blueprint Table I). Changelog

### `nf_update_derived_pub.m`
```matlab
function P = nf_update_derived_pub(P)
```
NF_UPDATE_DERIVED_PUB  Public wrapper — recompute aperture D and Rayleigh distance r_RD after changing P.M.  Also resets u bounds.

### `nf_usw_steer.m`
```matlab
function a = nf_usw_steer(theta, r, P)
```
NF_USW_STEER  Exact uniform spherical-wave (USW) near-field steering vector. a = nf_usw_steer(theta, r, P)

### `plot_fig.m`
```matlab
function plot_fig(x_vec, NMSE_db, NMSE_std, methods, title_str, xlabel_str, fig_name)
```
PLOT_FIG  Standard errorbar plot used by Figs 2, 3, 4, 7. plot_fig(x_vec, NMSE_db, NMSE_std, methods, title_str, xlabel_str, fig_name)

### `pregen.m`
```matlab
function [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR_dB, use_exact)
```
PREGEN  Pre-generate all Monte Carlo trial data before the estimation loop. Pre-generating data into cell arrays has two benefits:

### `save_results.m`
```matlab
function save_results(tag, NMSE_db, NMSE_std, x_vec, methods, varargin)
```
SAVE_RESULTS  Save simulation results to a .mat file for reproducibility. save_results(tag, NMSE_db, NMSE_std, x_vec, methods)

### `save_results_csv.m`
```matlab
function save_results_csv(csv_path, fig_id, fig_name, experiment_type, ...
```
SAVE_RESULTS_CSV  Append one row per method to the master results CSV. Called after each sweep point -- partial runs are preserved.

## Estimators (`estimators/`)

### `nf_bfsomp.m`
```matlab
function [theta_hat, r_hat, info] = nf_bfsomp(X, W, P)
```
NF_BFSOMP  Beam-Focused Simultaneous OMP (BF-SOMP). Implements Hussain, Abdallah & Eltawil,

### `nf_clkl.m`
```matlab
function [theta_hat, r_hat, p_hat, u_hat, info] = nf_clkl(R_hat, W, P)
```
NF_CLKL  Curvature-Learned KL Covariance Fitting estimator (proposed). IMPROVEMENTS in v5 (from paper audits):

### `nf_dfrft_nomp.m`
```matlab
function [theta_hat, r_hat, info] = nf_dfrft_nomp(X, ~, P)
```
NF_DFRFT_NOMP  DFrFT-based NOMP gridless hybrid-field channel estimator. Implements Yang et al. (IEEE Wireless Commun. Lett. 2024), Algorithm 1.

### `nf_music_tri.m`
```matlab
function [theta_hat, r_hat, info] = nf_music_tri(X, W_or_P, P_arg)
```
NF_MUSIC_TRI  Subarray MUSIC + geometric triangulation. Calling conventions (both supported):

### `nf_psomp.m`
```matlab
function [theta_hat, r_hat, info] = nf_psomp(R_hat, W, P)
```
NF_PSOMP  Polar-grid Simultaneous OMP (P-SOMP) -- v3. STEP 2 UPGRADE: Beam-depth range sampling (Hussain et al. TWC 2025).

### `nf_zhang.m`
```matlab
function [theta_hat, r_hat, info] = nf_zhang(X, ~, P)
```
NF_ZHANG  Dictionary-Learning OMP (DL-OMP) near-field channel estimator. Implements Zhang, Zhang & Eldar (IEEE Trans. Commun. 2024), Algorithm 1

## Figure Scripts (`figures/`)

### `run_fig2_nmse_snr.m`
```matlab
function run_fig2_nmse_snr(fast, use_par)
```
use_par : logical  (default false) -- pass true to use parfor across MC trials RUN_FIG2_NMSE_SNR  Fig. 2 -- NMSE vs SNR (primary figure).

### `run_fig3_nmse_nrf.m`
```matlab
function run_fig3_nmse_nrf(fast, use_par)
```
use_par : logical  (default false) -- pass true to use parfor across MC trials RUN_FIG3_NMSE_NRF  Fig. 3 -- NMSE vs N_RF (hybrid compression sweep).

### `run_fig4_nmse_n.m`
```matlab
function run_fig4_nmse_n(fast, use_par)
```
use_par : logical  (default false) — pass true to use parfor across MC trials RUN_FIG4_NMSE_N  Fig. 4 — NMSE vs N (snapshot sweep).

### `run_fig5_nearfar.m`
```matlab
function run_fig5_nearfar(fast, use_par)
```
use_par : logical (default false) -- enables parfor across MC trials RUN_FIG5_NEARFAR  Fig. 5 -- Near->Far transition sweep (Step 1 updated).

### `run_fig6_runtime.m`
```matlab
function run_fig6_runtime(fast, use_par)  %#ok<INUSD>
```
use_par is accepted but ignored -- Fig.6 always runs serially. Timing measurements must reflect single-core cost; parfor would

### `run_fig7_robustness.m`
```matlab
function run_fig7_robustness(fast, use_par)
```
use_par : logical (default false) -- enables parfor across MC trials RUN_FIG7_ROBUSTNESS  Fig. 7 -- Robustness: exact USW data, Fresnel estimators.

### `run_fig7b_source_robustness.m`
```matlab
function run_fig7b_source_robustness(fast, use_par)
```
RUN_FIG7B_SOURCE_ROBUSTNESS  Fig. 7b -- Source model robustness. New figure added in Step 7 of the TWC submission improvement plan.

### `run_fig8_vard.m`
```matlab
function run_fig8_vard(fast, use_par)
```
RUN_FIG8_VARD  Fig. 8 -- NMSE vs number of paths d. New figure added in Step 5 of the TWC submission improvement plan.

### `run_fig9_convergence.m`
```matlab
function run_fig9_convergence(fast, use_par)
```
RUN_FIG9_CONVERGENCE  Fig. 9 -- CL-KL convergence diagnostic. New figure added in Step 6 of the TWC submission improvement plan.

## Metrics & CRB (`metrics/`)

### `nf_crb.m`
```matlab
function [rmse_theta_crb, rmse_r_crb] = nf_crb(theta_true, r_true, p_true, N0, W, P)
```
NF_CRB  Stochastic Cramer-Rao Bound for compressed near-field observations. Computes the CRB for angle (theta) and range (r) estimation from the

### `nf_metrics.m`
```matlab
function [nmse, rmse_theta, rmse_r, fail_rate] = nf_metrics(...
```
NF_METRICS  Compute channel NMSE, angle/range RMSE, and failure rate. BUG-FIX (Bug 3 — matching cost robustness):

## Entry Points (root)

### `run_ablation.m`
```matlab
function run_ablation()
```
RUN_ABLATION  Ablation study for CL-KL (P2 response to reviewer Q1). Quantifies the individual contribution of each CL-KL design choice:

### `run_all.m`
```matlab
function run_all(mode)
```
RUN_ALL  Master runner -- executes all ten figures (Steps 1-7 of TWC plan). USAGE

### `setup.m`
*Script (no output arguments)*
SETUP  Add the simulation directory AND all subfolders to MATLAB path.

### `test_fixes.m`
*Script (no output arguments)*
TEST_FIXES  Bug-fix & algorithm-improvement verification suite.

### `test_modules.m`
*Script (no output arguments)*
TEST_MODULES  Smoke tests for all estimator functions.

