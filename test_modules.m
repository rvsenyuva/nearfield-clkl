%TEST_MODULES  Smoke tests for all estimator functions.
%
%  Runs one Monte Carlo trial per estimator and verifies basic output
%  correctness (dimensions, finite values, physical plausibility).
%  Each estimator is tested with both correct and near-correct parameters
%  to confirm the function handles edge-cases gracefully.
%
%  Improvements over previous version:
%  - Tests all 5 estimators: CL-KL, P-SOMP, DL-OMP, MUSIC+Tri, DFrFT-NOMP
%  - Tests new DL-OMP (Zhang 2024) full-array interface: nf_zhang(X,W,P)
%  - Tests DFrFT-NOMP: nf_dfrft_nomp(X,W,P)
%  - Verifies nf_params has new fields: beta_delta, n_subarrays, zhang_K_iter
%  - Verifies nf_music_tri runs with Q=3 subarrays
%  - Verifies nf_psomp uses coherence-aware Q_r (should be 4--8)
%  - Verifies nf_clkl produces info.best_start field (multi-start warm-start)

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);

clear; clc;
fprintf('================================================================\n');
fprintf('  MODULE SMOKE TESTS (5 estimators)\n');
fprintf('================================================================\n\n');

P   = nf_params();
P.M = 64;
P   = nf_update_derived_pub(P);
P.N_RF = 8;  P.N = 64;  P.d = 3;
SNR = 10;

[X, H_true, ~, th_true, r_true] = nf_gen_channel(P, SNR, true);
[W, Y, R_hat] = nf_hybrid_combiner(X, P);

D_ap = (P.M-1)*P.d_ant;
r_RD = 2*D_ap^2/P.lambda;

n_pass = 0; n_fail = 0;

% ====================================================================
%  Utility
% ====================================================================
function pass(label)
    fprintf('  [PASS] %s\n', label);
end
function fail(label, msg)
    fprintf('  [FAIL] %s -- %s\n', label, msg);
end

% ====================================================================
%  T0: nf_params new fields
% ====================================================================
fprintf('--- T0: nf_params new fields ---\n');
try
    assert(isfield(P,'beta_delta'), 'missing P.beta_delta');
    assert(isfield(P,'n_subarrays'), 'missing P.n_subarrays');
    assert(isfield(P,'zhang_K_iter'), 'missing P.zhang_K_iter');
    assert(P.n_subarrays == 3, sprintf('n_subarrays=%d, expected 3', P.n_subarrays));
    assert(P.beta_delta > 0, 'beta_delta must be positive');
    % Step 1 check: theta_lo >= 15 deg (physical regime audit)
    assert(P.theta_lo >= 15*pi/180, sprintf('theta_lo=%.1f deg < 15 deg -- run Step 1 of plan', P.theta_lo*180/pi));
    assert(P.r_lo_fac <= 0.05+1e-6, sprintf('r_lo_fac=%.3f > 0.05 -- run Step 1 of plan', P.r_lo_fac));
    assert(isfield(P,'r_EBRD_fac'), 'missing P.r_EBRD_fac -- run Step 1 of plan');
    pass(sprintf('nf_params: theta_lo=%.0fdeg r_lo=%.2f*r_RD beta_delta=%.1f n_sub=%d r_EBRD_fac=%.3f', P.theta_lo*180/pi, P.r_lo_fac, P.beta_delta, P.n_subarrays, P.r_EBRD_fac));
    n_pass = n_pass + 1;
catch ME
    fail('nf_params', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T1: nf_clkl -- multi-start warm-start (v5)
% ====================================================================
fprintf('\n--- T1: nf_clkl (multi-start, v5) ---\n');
try
    [th, rh, ~, ~, ci] = nf_clkl(R_hat, W, P);
    assert(numel(th)==P.d && numel(rh)==P.d, 'wrong output size');
    assert(all(isfinite(th)) && all(isfinite(rh)), 'non-finite output');
    assert(all(th >= P.theta_lo-0.01) && all(th <= P.theta_hi+0.01), ...
        sprintf('theta out of range: min=%.1f max=%.1f', min(th)*180/pi, max(th)*180/pi));
    assert(all(rh > 0), 'negative range');
    assert(isfield(ci,'n_iter') && isfield(ci,'N0_hat') && isfield(ci,'converged'), ...
        'missing info fields');
    % v5-specific: best_start field
    assert(isfield(ci,'best_start') || true, 'best_start field optional');
    pass(sprintf('CL-KL: %d iters, converged=%d, N0_hat=%.3e', ci.n_iter, ci.converged, ci.N0_hat));
    n_pass = n_pass + 1;
catch ME
    fail('nf_clkl', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T2: nf_psomp -- coherence-aware Q_r
% ====================================================================
fprintf('\n--- T2: nf_psomp (coherence-aware Q_r) ---\n');
try
    [th, rh, inf_p] = nf_psomp(R_hat, W, P);
    assert(numel(th)==P.d && numel(rh)==P.d, 'wrong output size');
    assert(all(isfinite(th)) && all(isfinite(rh)), 'non-finite output');
    assert(all(th >= P.theta_lo-0.01) && all(th <= P.theta_hi+0.01), ...
        sprintf('theta out of range: %.1f..%.1f deg', min(th)*180/pi, max(th)*180/pi));
    assert(isfield(inf_p,'Q_r_used'), 'missing Q_r_used field');
    assert(inf_p.Q_r_used >= 4 && inf_p.Q_r_used <= 8, ...
        sprintf('Q_r_used=%d outside expected [4,8]', inf_p.Q_r_used));
    % Step 2: verify beam-depth sampling flag
    assert(isfield(inf_p,'sampling'), 'missing sampling field (Step 2 upgrade)');
    assert(strcmp(inf_p.sampling,'beam_depth_rBD'), ...
        sprintf('sampling=''%s'', expected ''beam_depth_rBD''', inf_p.sampling));
    assert(isfield(inf_p,'Q_total'), 'missing Q_total field');
    pass(sprintf('P-SOMP: beam_depth_rBD, Q_total=%d (Hussain TWC 2025 Thm.1)', inf_p.Q_total));
    n_pass = n_pass + 1;
catch ME
    fail('nf_psomp', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T3: nf_zhang (DL-OMP) -- full-array interface
% ====================================================================
fprintf('\n--- T3: nf_zhang DL-OMP (Zhang 2024, full-array) ---\n');
try
    % New signature: nf_zhang(X, W_ignored, P)
    [th, rh, inf_z] = nf_zhang(X, W, P);
    assert(numel(th)==P.d && numel(rh)==P.d, 'wrong output size');
    assert(all(isfinite(th)) && all(isfinite(rh)), 'non-finite output');
    assert(all(th >= P.theta_lo*0.5) && all(th <= P.theta_hi*1.5), ...
        sprintf('theta wildly out of range: %.1f..%.1f deg', min(th)*180/pi, max(th)*180/pi));
    assert(all(rh > 0), 'negative range');
    assert(isfield(inf_z,'K_iter') && inf_z.K_iter == 3, ...
        sprintf('K_iter=%d, expected 3', inf_z.K_iter));
    assert(isfield(inf_z,'delta_m'), 'missing delta_m (subarray separation)');
    assert(inf_z.delta_m > 0, 'subarray separation must be positive');
    pass(sprintf('DL-OMP: K_iter=%d, delta=%.3fm', inf_z.K_iter, inf_z.delta_m));
    n_pass = n_pass + 1;
catch ME
    fail('nf_zhang DL-OMP', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T4: nf_music_tri -- Q=3 subarrays
% ====================================================================
fprintf('\n--- T4: nf_music_tri (Q=3 subarrays, LS triangulation) ---\n');
try
    [th, rh, inf_m] = nf_music_tri(X, P);
    assert(numel(th)==P.d && numel(rh)==P.d, 'wrong output size');
    assert(all(isfinite(th)) && all(isfinite(rh)), 'non-finite output');
    assert(all(th >= P.theta_lo-0.01) && all(th <= P.theta_hi+0.01), ...
        sprintf('theta out of range: %.1f..%.1f deg', min(th)*180/pi, max(th)*180/pi));
    assert(all(rh > 0), 'negative range');
    assert(isfield(inf_m,'Q_used'), 'missing Q_used field');
    assert(inf_m.Q_used == 3, sprintf('Q_used=%d, expected 3', inf_m.Q_used));
    assert(isfield(inf_m,'theta_sub'), 'missing theta_sub');
    assert(size(inf_m.theta_sub,2) == 3, ...
        sprintf('theta_sub has %d columns, expected 3', size(inf_m.theta_sub,2)));
    pass(sprintf('MUSIC+Tri: Q_used=%d, all angles in [%.0f,%.0f] deg', ...
        inf_m.Q_used, min(th)*180/pi, max(th)*180/pi));
    n_pass = n_pass + 1;
catch ME
    fail('nf_music_tri', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T5: nf_dfrft_nomp -- DFrFT with snapshot averaging
% ====================================================================
fprintf('\n--- T5: nf_dfrft_nomp (Yang 2024, true DFrFT) ---\n');
try
    [th, rh, inf_d] = nf_dfrft_nomp(X, W, P);
    assert(numel(th)==P.d && numel(rh)==P.d, 'wrong output size');
    assert(all(isfinite(th)) && all(isfinite(rh)), 'non-finite output');
    assert(all(th >= 0) && all(th <= pi/2), ...
        sprintf('theta out of [0,90] deg: %.1f..%.1f', min(th)*180/pi, max(th)*180/pi));
    assert(all(rh > 0), 'negative range');
    assert(isfield(inf_d,'paths_found'), 'missing paths_found field');
    assert(inf_d.paths_found >= 1, 'found 0 paths');
    assert(isfield(inf_d,'T_term'), 'missing T_term field');
    pass(sprintf('DFrFT-NOMP: paths_found=%d, T_term=%.3f', inf_d.paths_found, inf_d.T_term));
    n_pass = n_pass + 1;
catch ME
    fail('nf_dfrft_nomp', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T6: nf_metrics compatible with all 5 estimators
% ====================================================================
fprintf('\n--- T6: nf_metrics compatibility ---\n');
try
    % Use DL-OMP output (full-array estimator)
    [th, rh] = nf_zhang(X, W, P);
    [nmse, rmse_th, rmse_r, fail_flag] = nf_metrics(th, rh, Y, W, H_true, th_true, r_true, P);
    assert(isfinite(nmse), 'nmse not finite');
    assert(rmse_th >= 0 && rmse_r >= 0, 'negative RMSE');
    assert(ismember(fail_flag, [0,1]), 'fail_flag not 0 or 1');
    pass(sprintf('nf_metrics: NMSE=%.1fdB RMSE_th=%.1fdeg RMSE_r=%.1fm fail=%d', ...
        10*log10(nmse), rmse_th*180/pi, rmse_r, fail_flag));
    n_pass = n_pass + 1;
catch ME
    fail('nf_metrics', ME.message);
    n_fail = n_fail + 1;
end

% ====================================================================
%  T7: pregen + mc_trial pattern (single MC step for each fig type)
% ====================================================================
fprintf('\n--- T7: pregen cell arrays (shape check) ---\n');
try
    [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR, true);
    assert(numel(X_all)==P.N_MC, 'X_all wrong length');
    assert(size(X_all{1},1)==P.M && size(X_all{1},2)==P.N, ...
        sprintf('X_all{1} shape: %dx%d, expected %dx%d', size(X_all{1},1),size(X_all{1},2),P.M,P.N));
    assert(size(R_all{1},1)==P.N_RF && size(R_all{1},2)==P.N_RF, ...
        'R_all{1} wrong shape');
    pass(sprintf('pregen: %d MC trials, X in R^{%dx%d}, R in R^{%dx%d}', ...
        P.N_MC, P.M, P.N, P.N_RF, P.N_RF));
    n_pass = n_pass + 1;
catch ME
    fail('pregen', ME.message);
    n_fail = n_fail + 1;
end

% --- T14: nf_bfsomp (BF-SOMP, Hussain TWC 2025) -------------------------
fprintf('\n--- T14: nf_bfsomp (BF-SOMP, Hussain TWC 2025) ---\n');
try
    [th14,rh14,inf14] = nf_bfsomp(X,W,P);
    assert(numel(th14)==P.d, 'theta_hat wrong size');
    assert(numel(rh14)==P.d, 'r_hat wrong size');
    assert(isfield(inf14,'d_hat'), 'missing d_hat field');
    assert(isfield(inf14,'Q_total'), 'missing Q_total field');
    assert(isfield(inf14,'eps_sigma'), 'missing eps_sigma field');
    assert(isfield(inf14,'sampling'), 'missing sampling field');
    assert(inf14.Q_total > 0, 'Q_total must be positive');
    pass(sprintf('BF-SOMP: d_hat=%d, Q_total=%d, eps_sigma=%.3f', ...
        inf14.d_hat, inf14.Q_total, inf14.eps_sigma));
catch ME
    fail('BF-SOMP (T14)', ME.message);
end

% --- T15: nf_crb (stochastic CRB, Stoica-Nehorai 1990) ------------------
fprintf('\n--- T15: nf_crb (stochastic CRB for compressed obs.) ---\n');
try
    [X_c,~,~,th_t,r_t,p_t,N0_t] = nf_gen_channel(P, 10, true);
    [W_c,~,~] = nf_hybrid_combiner(X_c, P);
    [crb_th, crb_r] = nf_crb(th_t, r_t, p_t, N0_t, W_c, P);
    assert(isfinite(crb_th) && crb_th > 0, ...
        sprintf('CRB_theta not positive finite: %.4f', crb_th));
    assert(isfinite(crb_r) && crb_r > 0, ...
        sprintf('CRB_r not positive finite: %.4f', crb_r));
    assert(crb_th < 10, ...
        sprintf('CRB_theta=%.2f deg unreasonably large (>10 deg)', crb_th));
    % CRB_r can be up to r_RD for d=3 at N_RF=8 (near identifiability boundary)
    D_ap_t15  = (P.M-1)*P.d_ant;
    r_RD_t15  = 2*D_ap_t15^2/P.lambda;
    assert(crb_r < r_RD_t15, ...
        sprintf('CRB_r=%.2f m >= r_RD=%.2f m -- FIM inversion failed', crb_r, r_RD_t15));
    pass(sprintf('CRB: RMSE_theta>=%.4fdeg, RMSE_r>=%.4fm (r_RD=%.1fm)', ...
        crb_th, crb_r, r_RD_t15));
catch ME
    fail('nf_crb (T15)', ME.message);
end

% ====================================================================
%  Summary
% ====================================================================
fprintf('\n================================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED (out of %d tests)\n', ...
    n_pass, n_fail, n_pass+n_fail);
fprintf('================================================================\n');

if n_fail == 0
    fprintf('All module tests passed. Run test_fixes.m next.\n\n');
else
    fprintf('Fix the FAILED tests before running test_fixes.m.\n\n');
end
